/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(zewdie,PairZewdie)

#else

#ifndef LMP_PAIR_ZEWDIE_H
#define LMP_PAIR_ZEWDIE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairZewdie : public Pair {
 public:
  PairZewdie(LAMMPS *lmp);
  virtual ~PairZewdie();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);

 protected:
  enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

  //double cut_global;
  double **cut;

  //double gamma,upsilon,mu;   // Gay-Berne parameters
  double **shape1;           // per-type radii in x, y and z
  //double **shape2;           // per-type radii in x, y and z SQUARED
  //double *lshape;            // precalculation based on the shape
  //double **well;             // well depth scaling along each axis ^ -1.0/mu
  //double **epsilon,**sigma;  // epsilon and sigma values for atom-type pairs

  //Define parameters for Zewdie model
  double **epsilon0;
  double **sigma0;
  double epsilon000, epsiloncc2, epsilon220, epsilon222, epsilon224;
  double sigma000, sigmacc2, sigma220, sigma222, sigma224;

  int **form;  //SPHERE_SPEHRE, ELLIPSE_ELLIPSE, SPHERE_ELLIPSE, ELLIPSE_SPEHRE
  //double **lj1,**lj2,**lj3,**lj4;
  double **offset;
  //int *setwell;
  class AtomVecEllipsoid *avec;

  void allocate();

  double zewdie_analytic(const int i, const int j, double quatI[4], 
		           double quatJ[4], double *r12,
                           const double epsilon, const double sigma, double *fforce, double *ttor, double *rtor);

  double zewdie_lj(const int i,const int j, 
				       double *quat2, double *r12,
                       const double epsilon, const double sigma, double *fforce,
                       double *rtor);


  double zewdie_analytic_old(const int i, const int j, double quatI[4], 
		           double quatP1[4],double *r12,
                           const double rsq, double *fforce, double *ttor,
                           double *rtor);

  // orig gayberne
  //double zewdie_analytic(const int i, const int j, double a1[3][3],
  //                         double a2[3][3], double b1[3][3], double b2[3][3],
  //                         double g1[3][3], double g2[3][3], double *r12,
  //                         const double rsq, double *fforce, double *ttor,
  //                         double *rtor);
  //double zewdie_lj(const int i, const int j, double a1[3][3],
  //                   double b1[3][3],double g1[3][3],double *r12,
  //                   const double rsq, double *fforce, double *ttor);
  //void compute_eta_torque(double m[3][3], double m2[3][3],
  //                        double *s, double ans[3][3]);
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair zewdie requires atom style ellipsoid

Self-explanatory.

E: Pair zewdie requires atoms with same type have same shape

Self-explanatory.

E: Pair zewdie epsilon a,b,c coeffs are not all set

Each atom type involved in pair_style zewdie must
have these 3 coefficients set at least once.

E: Bad matrix inversion in mldivide3

This error should not occur unless the matrix is badly formed.

*/
