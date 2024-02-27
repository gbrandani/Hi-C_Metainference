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

//+PLUMEDOC ISDB_COLVAR HIC2
/*
Calculates the Hi-C contact map of a polymer, as described in:
Carstens, Nilges and Habeck, "Bayesian inference of chromatin structure ensembles from population-averaged contact data" PNAS 117.14 (2020)

\f[
f_{ij} = scale s( d_{ij} - d_0 )
\f]

where

\f[
s(d) = ( 1 - (d/r0)^6 ) / ( 1 - (d/r0)^12 )
\f]

The derivative is:

\f[
df_{ij}/dd_{ij} = scale/r0 (12x^11(1-x^6)-6x^5(1-x^12))/(1-x^12)^2
\f]

where

\f[
x = (d-d0)/r0
\f]

This collective variable calculates the contacts for a set of couple of atoms using the above definition.

Replica-Averaged simulations can be performed using HIC2s, \ref ENSEMBLE, \ref STATS and \ref RESTRAINT .
\ref METAINFERENCE can be activated using DOSCORE and the other relevant keywords.

In the plumed input file, it can be defined with, for example:

HIC2 ...
R_0=1.0 # the scale r0 (the smaller r0, the sharper the contact)
D_0=0.5 # a shift distance in nm (the default plumed unit of length)
ATOMS1=20,21 REFERENCE1=0.90
ATOMS2=37,38 REFERENCE2=0.15
LABEL=hic
... HIC2

the scale parameter is set from the "metainference" plumed command:
https://www.plumed.org/doc-v2.6/user-doc/html/_m_e_t_a_i_n_f_e_r_e_n_c_e.html

*/
//+ENDPLUMEDOC


class HIC2 :
  public MetainferenceBase
{
private:
  double         Const;
  double         r0;
  double         d0;
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
  explicit HIC2(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
  void calculate() override;
  void update() override;
};

PLUMED_REGISTER_ACTION(HIC2,"HIC2")

void HIC2::registerKeywords( Keywords& keys ) {
  componentsAreNotOptional(keys);
  MetainferenceBase::registerKeywords(keys);
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the couple of atoms involved in each of the bonds for which you wish to calculate the HIC2. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one contact will be "
           "calculated for each ATOMS keyword you specify.");
  keys.reset_style("ATOMS","atoms");
  keys.add("compulsory","R_0","1.","Add the scale r0. ");
  keys.add("compulsory","D_0","1.","Add the scaling factor to take into account concentration and other effects. ");
  keys.add("numbered","REFERENCE","Add an experimental value for each contact (useful for \\ref STATS).");
  keys.addOutputComponent("hic","default","the calculated # HIC2");
  keys.addOutputComponent("exp","REFERENCE","the experimental # HIC2");
}

HIC2::HIC2(const ActionOptions&ao):
  PLUMED_METAINF_INIT(ao),
  Const(1.),
  r0(1.),
  d0(1.),
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

  log.printf("  This is a fast implementation of HIC2.\n");

  // Read in R_0 constant
  parse("R_0", r0);
  if(r0<=0.) error("R_0 cannot be <= 0");

  // Read in shift factor d0
  parse("D_0", d0);
  if(d0<=0.) error("D_0 cannot be <= 0");

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
  log.printf("  R0 is %f. Scaling factor is %f.\n",r0,d0);
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
    error("  Cannot precompute HIC2 pointers, ndata!=getNumberOfAtoms()/2");
  }
  // GBB
  pntr_to_components.resize( ndata );
  for(unsigned i=0; i<ndata; i++) {
    std::string num; Tools::convert(i,num);
    pntr_to_components[i]=getPntrToComponent("hic-"+num);;
  }

}

void HIC2::calculate()
{

  const unsigned N=getNumberOfAtoms();
  vector<Vector> dHIC2(N/2, Vector{0.,0.,0.});

  /* HIC2 Calculations and forces */
  #pragma omp parallel num_threads(OpenMP::getNumThreads())
  {
    #pragma omp for
    for(unsigned r=0; r<N; r+=2)
    {
      const unsigned index=r/2;
      Vector  distance;
      if(pbc) distance    = pbcDistance(getPosition(r),getPosition(r+1));
      else    distance    = delta(getPosition(r),getPosition(r+1));
      const double dij    = sqrt(distance.modulo2());
      const double d      = dij-d0;
      const double ind    = 1.0/dij; // inverse distance
      const double x      = d/r0;
      const double x2     = x*x;
      const double x3     = x*x2;
      const double x5     = x2*x3;
      const double x6     = x3*x3;
      const double x11    = x5*x6;
      const double x12    = x6*x6;
      const double one_x6     = 1.0 - x6;
      const double one_x12    = 1.0 - x12;
      const double inv_one_x12= 1.0/one_x12;
      const double epsilon= 0.000001;

      double hic;
      if( x>1.0-epsilon && x<1.0+epsilon) hic = 0.5;
      else                                hic = one_x6*inv_one_x12;

      // derivative
      double der_term;
      if( x>1.0-epsilon && x<1.0+epsilon) der_term = -1.5;
      else                                der_term = (12.0*x11*one_x6-6.0*x5*one_x12)*inv_one_x12*inv_one_x12;
      dHIC2[index] = - distance * ind * der_term / r0;

      string num; Tools::convert(index,num);
      // this is one step that makes metainference scale quadratically
      // can we precalculate the pointers to the values at the beginning and store them into an array?
      //Value* val=getPntrToComponent("hic-"+num);
      //val->set(hic);
      //if(!getDoScore()) {
      //  setBoxDerivatives(val, Tensor(distance,dHIC2[index]));
      //  setAtomsDerivatives(val, r,  dHIC2[index]);
      //  setAtomsDerivatives(val, r+1, -dHIC2[index]);
      //} else setCalcData(index, hic);
      pntr_to_components[index]->set(hic);
      if(!getDoScore()) {
        setBoxDerivatives(pntr_to_components[index], Tensor(distance,dHIC2[index]));
        // GBB derivatives of a value now store only the two corresponding atoms, not all of them
        //setAtomsDerivatives(pntr_to_components[index], r,  dHIC2[index]);
        //setAtomsDerivatives(pntr_to_components[index], r+1, -dHIC2[index]);
        setAtomsDerivatives(pntr_to_components[index], 0,  dHIC2[index]); // GBB
        setAtomsDerivatives(pntr_to_components[index], 1, -dHIC2[index]); // GBB
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
      const Vector der = dHIC2[index]*getMetaDer(index);
      dervir += Tensor(distance, der);
      setAtomsDerivatives(val, r,  der);
      setAtomsDerivatives(val, r+1, -der);
    }
    setBoxDerivatives(val, dervir);
  }
}

void HIC2::update() {
  // write status file
  if(getWstride()>0&& (getStep()%getWstride()==0 || getCPT()) ) writeStatus();
}

}
}
