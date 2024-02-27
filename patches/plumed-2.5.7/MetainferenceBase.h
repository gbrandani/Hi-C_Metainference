/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#ifndef __PLUMED_isdb_MetainferenceBase_h
#define __PLUMED_isdb_MetainferenceBase_h

#include "core/ActionWithValue.h"
#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "core/PlumedMain.h"
#include "tools/Random.h"
#include "tools/OpenMP.h"

#define PLUMED_METAINF_INIT(ao) Action(ao),MetainferenceBase(ao)

namespace PLMD {
namespace isdb {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new ISDB Metainference actions, within it there is
information as to how to go about implementing a new Metainference action.
*/

class MetainferenceBase :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue
{
private:
  std::vector<double> forces;
  std::vector<double> forcesToApply;

  // activate metainference
  bool doscore_;
  unsigned write_stride_;
  // number of experimental data
  unsigned narg;
  // experimental data
  std::vector<double> parameters;
  // metainference derivatives
  std::vector<double> metader_;
  // vector of back-calculated experimental data
  std::vector<double> calc_data_;

  // noise type
  unsigned noise_type_;
  enum { GAUSS, MGAUSS, OUTLIERS, MOUTLIERS, GENERIC };
  unsigned gen_likelihood_;
  enum { LIKE_GAUSS, LIKE_LOGN };
  bool   doscale_;
  unsigned scale_prior_;
  enum { SC_GAUSS, SC_FLAT };
  double scale_;
  double scale_mu_;
  double scale_min_;
  double scale_max_;
  double Dscale_;
  // scale is data scaling factor
  // noise type
  unsigned offset_prior_;
  bool   dooffset_;
  double offset_;
  double offset_mu_;
  double offset_min_;
  double offset_max_;
  double Doffset_;
  // scale and offset regression
  bool doregres_zero_;
  int  nregres_zero_;
  // sigma is data uncertainty
  std::vector<double> sigma_;
  std::vector<double> sigma_min_;
  std::vector<double> sigma_max_;
  std::vector<double> Dsigma_;
  // sigma_mean is uncertainty in the mean estimate
  std::vector<double> sigma_mean2_;
  // this is the estimator of the mean value per replica for generic metainference
  std::vector<double> ftilde_;
  double Dftilde_;

  // temperature in kbt
  double   kbt_;

  // Monte Carlo stuff
  std::vector<Random> random;
  unsigned MCsteps_;
  unsigned MCstride_;
  long unsigned MCaccept_;
  long unsigned MCacceptScale_;
  long unsigned MCacceptFT_;
  long unsigned MCtrial_;
  unsigned MCchunksize_;

  // output
  Value*   valueScore;
  Value*   valueScale;
  Value*   valueOffset;
  Value*   valueAccept;
  Value*   valueAcceptScale;
  Value*   valueAcceptFT;
  std::vector<Value*> valueSigma;
  std::vector<Value*> valueSigmaMean;
  std::vector<Value*> valueFtilde;

  // restart
  std::string status_file_name_;
  OFile    sfile_;

  // others
  bool     firstTime;
  std::vector<bool> firstTimeW;
  bool     master;
  bool     do_reweight_;
  unsigned do_optsigmamean_;
  unsigned nrep_;
  unsigned replica_;

  // selector
  unsigned nsel_;
  std::string selector_;
  unsigned iselect;

  // optimize sigma mean
  std::vector< std::vector < std::vector <double> > > sigma_mean2_last_;
  unsigned optsigmamean_stride_;

  // average weights
  double decay_w_;
  std::vector< std::vector <double> >  average_weights_;

  double getEnergyMIGEN(const std::vector<double> &mean, const std::vector<double> &ftilde, const std::vector<double> &sigma,
                        const double scale, const double offset);
  double getEnergySP(const std::vector<double> &mean, const std::vector<double> &sigma,
                     const double scale, const double offset);
  double getEnergySPE(const std::vector<double> &mean, const std::vector<double> &sigma,
                      const double scale, const double offset);
  double getEnergyGJ(const std::vector<double> &mean, const std::vector<double> &sigma,
                     const double scale, const double offset);
  double getEnergyGJE(const std::vector<double> &mean, const std::vector<double> &sigma,
                      const double scale, const double offset);
  void   setMetaDer(const unsigned index, const double der);
  double getEnergyForceSP(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  double getEnergyForceSPE(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  double getEnergyForceGJ(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  double getEnergyForceGJE(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  double getEnergyForceMIGEN(const std::vector<double> &mean, const std::vector<double> &dmean_x, const std::vector<double> &dmean_b);
  double getCalcData(const unsigned index);
  void get_weights(double &fact, double &var_fact);
  void replica_averaging(const double fact, std::vector<double> &mean, std::vector<double> &dmean_b);
  void get_sigma_mean(const double fact, const double var_fact, const std::vector<double> &mean);
  void do_regression_zero(const std::vector<double> &mean);
  void doMonteCarlo(const std::vector<double> &mean);


public:
  static void registerKeywords( Keywords& keys );
  explicit MetainferenceBase(const ActionOptions&);
  ~MetainferenceBase();
  void Initialise(const unsigned input);
  void Selector();
  unsigned getNarg();
  void setNarg(const unsigned input);
  void setParameters(const std::vector<double>& input);
  void setParameter(const double input);
  void setCalcData(const unsigned index, const double datum);
  void setCalcData(const std::vector<double>& data);
  bool getDoScore();
  unsigned getWstride();
  double getScore();
  void setScore(const double score);
  void setDerivatives();
  double getMetaDer(const unsigned index);
  void writeStatus();
  void turnOnDerivatives();
  unsigned getNumberOfDerivatives();
  void lockRequests();
  void unlockRequests();
  void calculateNumericalDerivatives( ActionWithValue* a ); // GBB FIXME this cannot be used
  void apply();
  void setArgDerivatives(Value *v, const double &d);
  void addForcesToApply(std::vector<double>& forcesToApply, std::vector<double>& forces, const unsigned i_component); // GBB
  void setAtomsDerivatives(Value*v, const unsigned i, const Vector&d);
  void setBoxDerivatives(Value*v, const Tensor&d);
};

inline
void MetainferenceBase::setNarg(const unsigned input)
{
  narg = input;
}

inline
bool MetainferenceBase::getDoScore()
{
  return doscore_;
}

inline
unsigned MetainferenceBase::getWstride()
{
  return write_stride_;
}

inline
unsigned MetainferenceBase::getNarg()
{
  return narg;
}

inline
void MetainferenceBase::setMetaDer(const unsigned index, const double der)
{
  metader_[index] = der;
}

inline
double MetainferenceBase::getMetaDer(const unsigned index)
{
  return metader_[index];
}

inline
double MetainferenceBase::getCalcData(const unsigned index)
{
  return calc_data_[index];
}

inline
void MetainferenceBase::setCalcData(const unsigned index, const double datum)
{
  calc_data_[index] = datum;
}

inline
void MetainferenceBase::setCalcData(const std::vector<double>& data)
{
  for(unsigned i=0; i<data.size(); i++) calc_data_[i] = data[i];
}

inline
void MetainferenceBase::setParameters(const std::vector<double>& input) {
  for(unsigned i=0; i<input.size(); i++) parameters.push_back(input[i]);
}

inline
void MetainferenceBase::setParameter(const double input) {
  parameters.push_back(input);
}

inline
void MetainferenceBase::setScore(const double score) {
  valueScore->set(score);
}

inline
void MetainferenceBase::setDerivatives() {
  // GBB in the original version, the high number of (unused) derivatives makes the code slow
  // Get appropriate number of derivatives
  // Derivatives are first for arguments and then for atoms
  unsigned nder;
  unsigned nder_per_value;
  // GBB FIXME this is ok for HIC, but may need modifications for the other collective variables
  const unsigned natoms_per_value = 2;
  // GBB derivatives should correspond to only 2 atoms
  if( getNumberOfAtoms()>0 ) {
    nder           = 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
    nder_per_value = 3*natoms_per_value + 9 + getNumberOfArguments();
    //nder_per_value = 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
  } else {
    nder           = getNumberOfArguments();
    nder_per_value = getNumberOfArguments();
  }

  // Resize all derivative arrays
  // GBB
  //forces.resize( nder );
  forces.resize( nder_per_value ); // this should correspond just to 2 atoms
  forcesToApply.resize( nder ); // this sums over all the derivatives
  //for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->resizeDerivatives(nder);
  for(int i=0; i<getNumberOfComponents(); ++i) getPntrToComponent(i)->resizeDerivatives(nder_per_value);
}

inline
void MetainferenceBase::turnOnDerivatives() {
  ActionWithValue::turnOnDerivatives();
}

inline
unsigned MetainferenceBase::getNumberOfDerivatives() {
  if( getNumberOfAtoms()>0 ) {
    // GBB
    // FIXME this is ok for HIC, but may need modifications for the other collective variables
    //return 3*getNumberOfAtoms() + 9 + getNumberOfArguments();
    const unsigned nat=2;
    return 3*nat + 9 + getNumberOfArguments();
  }
  return getNumberOfArguments();
}

inline
void MetainferenceBase::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

inline
void MetainferenceBase::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

// GBB FIXME cannot use this
inline
void MetainferenceBase::calculateNumericalDerivatives( ActionWithValue* a=NULL ) {
  error("  MetainferenceBase::calculateNumericalDerivative() is broken, fix before using\\"); // GBB
  if( getNumberOfArguments()>0 ) {
    ActionWithArguments::calculateNumericalDerivatives( a );
  }
  if( getNumberOfAtoms()>0 ) {
    Matrix<double> save_derivatives( getNumberOfComponents(), getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) if(getPntrToComponent(j)->hasDerivatives()) save_derivatives(j,i)=getPntrToComponent(j)->getDerivative(i);
    }
    calculateAtomicNumericalDerivatives( a, getNumberOfArguments() );
    for(int j=0; j<getNumberOfComponents(); ++j) {
      for(unsigned i=0; i<getNumberOfArguments(); ++i) if(getPntrToComponent(j)->hasDerivatives()) getPntrToComponent(j)->addDerivative( i, save_derivatives(j,i) );
    }
  }
}

inline
void MetainferenceBase::apply() {
  // GBB this step is very slow when there are many components
  bool wasforced=false; forcesToApply.assign(forcesToApply.size(),0.0);
  for(int i=0; i<getNumberOfComponents(); ++i) {
    // GBB this was initially very slow
    // here you take the derivatives stored in the Value and pass them to forces
    //if(getPntrToComponent(i)->hasDerivatives()) {
    //}
    if( getPntrToComponent(i)->applyForce( forces ) ) {
      wasforced=true;
      // GBB very slow, too many forces in original version
      // in the new version forcesToApply stores all the forces to all the atoms (as before)
      // but forces only stores the forces to the 2 atoms of the corresponding Value
      //for(unsigned i=0; i<forces.size(); ++i) forcesToApply[i]+=forces[i];
      // get the indeces to which you want to apply the forces
      addForcesToApply(forcesToApply, forces, i); // GBB
    }
  }
  if( wasforced ) {
    addForcesOnArguments( forcesToApply );
    if( getNumberOfAtoms()>0 ) setForcesOnAtoms( forcesToApply, getNumberOfArguments() );
  }
}

inline
void MetainferenceBase::setArgDerivatives(Value *v, const double &d) {
  v->addDerivative(0,d);
}

// GBB
inline
void MetainferenceBase::addForcesToApply(std::vector<double>& forcesToApply, std::vector<double>& forces, const unsigned i_component) {
  const unsigned noa=getNumberOfArguments();
  const unsigned i_atom1 = i_component*2;
  const unsigned i_atom2 = i_component*2+1;
  forcesToApply[noa+3*i_atom1+0]+=forces[noa+3*0+0];
  forcesToApply[noa+3*i_atom1+1]+=forces[noa+3*0+1];
  forcesToApply[noa+3*i_atom1+2]+=forces[noa+3*0+2];
  forcesToApply[noa+3*i_atom2+0]+=forces[noa+3*1+0];
  forcesToApply[noa+3*i_atom2+1]+=forces[noa+3*1+1];
  forcesToApply[noa+3*i_atom2+2]+=forces[noa+3*1+2];
}

inline
void MetainferenceBase::setAtomsDerivatives(Value*v, const unsigned i, const Vector&d) {
  const unsigned noa=getNumberOfArguments();
  v->addDerivative(noa+3*i+0,d[0]);
  v->addDerivative(noa+3*i+1,d[1]);
  v->addDerivative(noa+3*i+2,d[2]);
}

inline
void MetainferenceBase::setBoxDerivatives(Value* v,const Tensor&d) {
  const unsigned noa=getNumberOfArguments();
  // GBB in the new version, each takes care of only 2 atoms
  // FIXME this is ok for HIC, but may need modifications for the other collective variables
  //const unsigned nat=getNumberOfAtoms();
  const unsigned nat=2;
  v->addDerivative(noa+3*nat+0,d(0,0));
  v->addDerivative(noa+3*nat+1,d(0,1));
  v->addDerivative(noa+3*nat+2,d(0,2));
  v->addDerivative(noa+3*nat+3,d(1,0));
  v->addDerivative(noa+3*nat+4,d(1,1));
  v->addDerivative(noa+3*nat+5,d(1,2));
  v->addDerivative(noa+3*nat+6,d(2,0));
  v->addDerivative(noa+3*nat+7,d(2,1));
  v->addDerivative(noa+3*nat+8,d(2,2));
}


}
}

#endif

