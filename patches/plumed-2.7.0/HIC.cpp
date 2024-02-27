/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "MetainferenceBase.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

#ifdef __PLUMED_HAS_GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_COLVAR HIC
/*
Calculates the Hi-C contact map of a polymer, as described in:
Carstens, Nilges and Habeck, "Bayesian inference of chromatin structure ensembles from population-averaged contact data" PNAS 117.14 (2020)

\f[
f_{ij} = scale s( d_c - d_{ij} )
\f]

where

\f[
s(d) = 0.5 ( 1 + gamma d / \sqrt{ 1 + gamma^2 d^2 } )
\f]

The derivative is:

\f[
df_{ij}/dd_{ij} = - 0.5 scale {\bf d_{ij}}/d_{ij} dg/dd(d_c-d_{ij})
\f]

where

\f[
dg/dd = gamma/\sqrt{1+gamma^2 d^2} - gamma^3 d^2 / (1+gamma^2 d^2)^{3/2}
\f]

This collective variable calculates the contacts for a set of couple of atoms using the above definition.

Replica-Averaged simulations can be performed using HICs, \ref ENSEMBLE, \ref STATS and \ref RESTRAINT .
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

In the plumed input file, it can be defined with, for example:

HIC ...
GAMMA=5.0 # steepness of the step function (the higher the sharper, in the limit of gamma=+inf you get a step function)
DC=1.0 # distance cutoff in nm (the default plumed unit of length)
ATOMS1=20,21 REFERENCE1=0.90
ATOMS2=37,38 REFERENCE2=0.15
LABEL=hic
... HIC

the scale parameter is set from the "metainference" plumed command:
https://www.plumed.org/doc-v2.6/user-doc/html/_m_e_t_a_i_n_f_e_r_e_n_c_e.html

*/
//+ENDPLUMEDOC


class HIC :
  public MetainferenceBase
{
private:
  double         Const;
  double         gamma;
  double         dc;
  vector<double> reference;
  bool           pbc;

  // GBB
  // store pointers to components
  std::vector<Value*> pntr_to_components;

#ifdef __PLUMED_HAS_GSL
/// Auxiliary class to delete a gsl_vector.
/// If used somewhere else we can move it.
  struct gsl_vector_deleter {
    void operator()(gsl_vector* p) {
      gsl_vector_free(p);
    }
  };

/// unique_ptr to a gsl_vector.
/// Gets deleted when going out of scope.
  typedef std::unique_ptr<gsl_vector,gsl_vector_deleter> gsl_vector_unique_ptr;

/// Auxiliary class to delete a gsl_matrix.
/// If used somewhere else we can move it.
  struct gsl_matrix_deleter {
    void operator()(gsl_matrix* p) {
      gsl_matrix_free(p);
    }
  };

/// unique_ptr to a gsl_matrix.
/// Gets deleted when going out of scope.
  typedef std::unique_ptr<gsl_matrix,gsl_matrix_deleter> gsl_matrix_unique_ptr;
#endif


public:
  explicit HIC(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(HIC,"HIC")

void HIC::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the couple of atoms involved in each of the bonds for which you wish to calculate the HIC. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one contact will be "
           "calculated for each ATOMS keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("compulsory","GAMMA","1.","Add the product of the gamma for the bond. ");
  keys.add("compulsory","DC","1.","Add the scaling factor to take into account concentration and other effects. ");
  keys.add("numbered","REFERENCE","Add an experimental value for each contact (useful for \\ref STATS).");
  keys.addOutputComponent("hic","default","the calculated # HIC");
  keys.addOutputComponent("exp","REFERENCE","the experimental # HIC");
}

HIC::HIC(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  Const(1.),
  gamma(1.),
  dc(1.),
  pbc(true)
{
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  // Read in the atoms
  vector<AtomNumber> t, atoms;
  for(int i=1;; ++i ) {
    parseAtomList("ATOMS", i, t );
    if( t.empty() ) break;
    if( t.size()!=2 ) {
      std::string ss; Tools::convert(i,ss);
      error("ATOMS" + ss + " keyword has the wrong number of atoms");
    }
    atoms.push_back(t[0]);
    atoms.push_back(t[1]);
    t.resize(0);
  }

  const unsigned ndata = atoms.size()/2;

  log.printf("  This is a fast implementation of HIC.\n");

  // Read in GAMMAAGNETIC constants
  parse("GAMMA", gamma);
  if(gamma<=0.) error("GAMMA cannot be <= 0");

  // Read in SCALING factors
  parse("DC", dc);
  if(dc<=0.) error("DC cannot be <= 0");

  // Optionally add an experimental value
  reference.resize( ndata );
  unsigned ntarget=0;
  for(unsigned i=0; i<ndata; ++i) {
    if( !parseNumbered( "REFERENCE", i+1, reference[i] ) ) break;
    ntarget++;
  }
  bool addexp=false;
  if(ntarget!=ndata && ntarget!=0) error("found wrong number of REFERENCE values");
  if(ntarget==ndata) addexp=true;
  if(getDoScore()&&!addexp) error("with DOSCORE you need to set the REFERENCE values");

  // Output details of all contacts
  log.printf("  Gamma is %f. Scaling factor is %f.\n",gamma,dc);
  for(unsigned i=0; i<ndata; ++i) {
    log.printf("  The %uth contact is calculated from atoms : %d %d.", i+1, atoms[2*i].serial(), atoms[2*i+1].serial());
    if(addexp) log.printf(" Experimental contact is %f.", reference[i]);
    log.printf("\n");
  }

  log<<"  Bibliography "
     <<plumed.cite("Camilloni C, Vendruscolo M, J. Phys. Chem. B 119, 653 (2015)")
     <<plumed.cite("Camilloni C, Vendruscolo M, Biochemistry 54, 7470 (2015)");
  log<< plumed.cite("Bonomi, Camilloni, Bioinformatics, 33, 3999 (2017)") << "\n";


  if(!getDoScore()) {
    for(unsigned i=0; i<ndata; i++) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("hic-"+num);
      componentIsNotPeriodic("hic-"+num);
    }
    if(addexp) {
      for(unsigned i=0; i<ndata; i++) {
        std::string num; Tools::convert(i,num);
        addComponent("exp-"+num);
        componentIsNotPeriodic("exp-"+num);
        Value* comp=getPntrToComponent("exp-"+num);
        comp->set(reference[i]);
      }
    }
  } else {
    for(unsigned i=0; i<ndata; i++) {
      std::string num; Tools::convert(i,num);
      addComponentWithDerivatives("hic-"+num);
      componentIsNotPeriodic("hic-"+num);
    }
    for(unsigned i=0; i<ndata; i++) {
      std::string num; Tools::convert(i,num);
      addComponent("exp-"+num);
      componentIsNotPeriodic("exp-"+num);
      Value* comp=getPntrToComponent("exp-"+num);
      comp->set(reference[i]);
    }
  }

  requestAtoms(atoms, false);
  if(getDoScore()) {
    setParameters(reference);
    Initialise(reference.size());
  }
  setDerivatives();
  checkRead();

  // precompute hic pointers
  if(ndata!=getNumberOfAtoms()/2) {
    error("  Cannot precompute HIC pointers, ndata!=getNumberOfAtoms()/2");
  }
  // GBB
  pntr_to_components.resize( ndata );
  for(unsigned i=0; i<ndata; i++) {
    std::string num; Tools::convert(i,num);
    pntr_to_components[i]=getPntrToComponent("hic-"+num);;
  }

}

void HIC::calculate()
{

  const unsigned N=getNumberOfAtoms();
  vector<Vector> dHIC(N/2, Vector{0.,0.,0.});

  /* HIC Calculations and forces */
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    #pragma omp for
    for(unsigned r=0; r<N; r+=2)
    {
      const unsigned index=r/2;
      Vector  distance;
      if(pbc) distance   = pbcDistance(getPosition(r),getPosition(r+1));
      else    distance   = delta(getPosition(r),getPosition(r+1));
      const double dij   = sqrt(distance.modulo2());
      const double d     = dc-dij;
      const double ind   = 1.0/dij; // inverse distance
      const double gamma_d       = gamma*d;
      const double inv_sqrt_term = 1.0 / sqrt(1.0+gamma_d*gamma_d);

      const double hic   = 0.5 * ( 1.0 + gamma_d * inv_sqrt_term );

      // derivative
      const double der_term = gamma*inv_sqrt_term - gamma*gamma_d*gamma_d*inv_sqrt_term*inv_sqrt_term*inv_sqrt_term;
      dHIC[index] = 0.5 * distance * ind * der_term;

      string num; Tools::convert(index,num);
      // this is one step that makes metainference scale quadratically
      // can we precalculate the pointers to the values at the beginning and store them into an array?
      //Value* val=getPntrToComponent("hic-"+num);
      //val->set(hic);
      //if(!getDoScore()) {
      //  setBoxDerivatives(val, Tensor(distance,dHIC[index]));
      //  setAtomsDerivatives(val, r,  dHIC[index]);
      //  setAtomsDerivatives(val, r+1, -dHIC[index]);
      //} else setCalcData(index, hic);
      pntr_to_components[index]->set(hic);
      if(!getDoScore()) {
        setBoxDerivatives(pntr_to_components[index], Tensor(distance,dHIC[index]));
        // GBB derivatives of a value now store only the two corresponding atoms, not all of them
        //setAtomsDerivatives(pntr_to_components[index], r,  dHIC[index]);
        //setAtomsDerivatives(pntr_to_components[index], r+1, -dHIC[index]);
        setAtomsDerivatives(pntr_to_components[index], 0,  dHIC[index]); // GBB
        setAtomsDerivatives(pntr_to_components[index], 1, -dHIC[index]); // GBB
      } else setCalcData(index, hic);
    }
  }

  if(getDoScore()) {
    /* Metainference */
    Tensor dervir;
    double score = getScore();
    setScore(score);

    /* calculate final derivatives */
    Value* val=getPntrToComponent("score");
    for(unsigned r=0; r<N; r+=2)
    {
      const unsigned index=r/2;
      Vector  distance;
      if(pbc) distance   = pbcDistance(getPosition(r),getPosition(r+1));
      else    distance   = delta(getPosition(r),getPosition(r+1));
      const Vector der = dHIC[index]*getMetaDer(index);
      dervir += Tensor(distance, der);
      setAtomsDerivatives(val, r,  der);
      setAtomsDerivatives(val, r+1, -der);
    }
    setBoxDerivatives(val, dervir);
  }
}

void HIC::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

}
}
