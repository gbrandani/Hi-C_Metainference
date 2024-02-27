/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Craig Tenney (University of Notre Dame)
     support for bond and angle restraints by Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_restrain_orient.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "respa.h"
#include "input.h"
#include "math_const.h"
#include "math_vector.h"
#include "memory.h"
#include "error.h"
#include "atom_vec_ellipsoid.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{ALIGN_F,ALIGN_V,ALIGN_U, THETA1, THETA2, PHI};

#define TOLERANCE 0.05
#define SMALL 0.001
#define DELTA 1

/* ---------------------------------------------------------------------- */

FixRestrainOrient::FixRestrainOrient(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix restrain command");

  scalar_flag = 1;
  vector_flag = 1;

  global_freq = 1;
  extscalar = 1;
  //respa_level_support = 1;
  ilevel_respa = 0;

  // parse args

  nrestrain = maxrestrain = 0;
  rstyle = NULL;
  rtype = NULL;
  ids = NULL;
  kstart = kstop = NULL;
  target = NULL;

  int iarg = 6;

  //loop through the memory
  while (iarg < narg) {
    if (nrestrain == maxrestrain) {
      maxrestrain += DELTA;
      memory->grow(rstyle,maxrestrain,"restrain:rstyle");
      memory->grow(ids,maxrestrain,4,"restrain:ids");
      memory->grow(rtype,maxrestrain,"restrain:rtype");
      memory->grow(kstart,maxrestrain,"restrain:kstart");
      memory->grow(kstop,maxrestrain,"restrain:kstop");
      memory->grow(target,maxrestrain,"restrain:target");
      memory->grow(angles,maxrestrain,"restrain:angles");
      //memory->grow(cos_target,maxrestrain,"restrain:cos_target");
      //memory->grow(sin_target,maxrestrain,"restrain:sin_target");
    }
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix restrain command");
      
      //check to make sure parameters line up
      if ((strcmp(arg[3],"f") == 0) || 
      (strcmp(arg[3],"v") == 0) ||
      (strcmp(arg[3],"u") == 0)) {
          if (strcmp(arg[3],"f") == 0)
          rstyle[nrestrain] = ALIGN_F;
          else if (strcmp(arg[3],"v") == 0) 
          rstyle[nrestrain] = ALIGN_V;
          else if (strcmp(arg[3],"u") == 0) 
          rstyle[nrestrain] = ALIGN_U;
          else error->all(FLERR,"Illegal fix restrain command");
          
          ids[nrestrain][0] = force->inumeric(FLERR,arg[4]);
          ids[nrestrain][1] = force->inumeric(FLERR,arg[5]);
      }
      else error->all(FLERR,"Illegal fix restrain command");

      //align the f,v,u vector of an atom to the bond vector (from angle_orient)
      //choose which angle you wish to restrain 
      if (strcmp(arg[iarg],"theta1") == 0) {
          rtype[nrestrain] = THETA1;
          target[nrestrain] = 180.0 - force->numeric(FLERR,arg[iarg+3]);
      }
      else if (strcmp(arg[iarg],"theta2") == 0) {
          rtype[nrestrain] = THETA2;
          target[nrestrain] = 180.0 - force->numeric(FLERR,arg[iarg+3]);
      }
      else if (strcmp(arg[iarg],"phi") == 0) {
          rtype[nrestrain] = PHI;
          target[nrestrain] = force->numeric(FLERR,arg[iarg+3]);
      }
      else error->all(FLERR,"Illegal fix restrain command");

      kstart[nrestrain] = force->numeric(FLERR,arg[iarg+1]);
      kstop[nrestrain] = force->numeric(FLERR,arg[iarg+2]);
      target[nrestrain] *= MY_PI/180.0;
      iarg += 4; 

    nrestrain++;
  }
  
  size_vector = nrestrain; //one element for each restraint, for the angle
  //angles = new double[size_vector];

  // require atom map to lookup atom IDs

  if (atom->map_style == 0)
    error->all(FLERR,"Fix restrain requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

FixRestrainOrient::~FixRestrainOrient()
{
  memory->destroy(rstyle);
  memory->destroy(ids);
  memory->destroy(rtype);
  memory->destroy(kstart);
  memory->destroy(kstop);
  memory->destroy(target);
  //memory->destroy(cos_target);
  //memory->destroy(sin_target);
  memory->destroy(angles);
}

/* ---------------------------------------------------------------------- */

int FixRestrainOrient::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRestrainOrient::init()
{
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    //if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixRestrainOrient::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixRestrainOrient::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRestrainOrient::post_force(int vflag)
{
  energy = 0.0;

  for (int m = 0; m < nrestrain; m++)
    if ((rstyle[m] == ALIGN_F) || (rstyle[m] == ALIGN_V) || (rstyle[m] == ALIGN_U))
      restrain_alignment(m);
}

/* ---------------------------------------------------------------------- */

void FixRestrainOrient::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRestrainOrient::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   constrain f,v,u of any two sites using the theta1, theta2, and phi 
   angles
---------------------------------------------------------------------- */

void FixRestrainOrient::restrain_alignment(int m)
{
  int i,j;
  double dx,dy,dz,fbond;
  //double rsq,r,dr,rk;

  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  i = atom->map(ids[m][0]);
  j = atom->map(ids[m][1]);
  if ((i == -1) || (j == -1))  {
      char str[128];
      sprintf(str,
              "restrain_orient atom %d %d missing on proc %d at step " BIGINT_FORMAT,
              ids[m][0],ids[m][1],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
  }
 

  /* ----------------------------------------------
   Allignment torque
   ---------------------------------------------- */
  double u0[3],uI[3],uJ[3],r[3], rn;
  double fterm1[3],fterm2[3],force1[3],force3[3];
  double *quatI, *quatJ;
  double quatItmp[4], quatJtmp[4];
  AtomVecEllipsoid *avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  int *ellipsoid = atom->ellipsoid;
  double U,theta1, theta2, phi,dtheta1, dtheta2, dphi;
  double inverse;
  double a0,a1,a2, dUda0, dUda1, dUda2;
  double dUdfI[3], dUdfJ[3];
  double torq[3], torqI[3], torqJ[3];
  
  if (rstyle[m] == ALIGN_F){
    u0[0] = 1; u0[1] = 0; u0[2] = 0; //f
  }
  else if (rstyle[m] == ALIGN_V){
    u0[0] = 0; u0[1] = 1; u0[2] = 0; //v
  }
  else if (rstyle[m] == ALIGN_U){
    u0[0] = 0; u0[1] = 0; u0[2] = 1; //u
  }

  //get rIJ
  dx = x[i][0] - x[j][0]; // get vector between i and j
  dy = x[i][1] - x[j][1];
  dz = x[i][2] - x[j][2];
  r[0] = dx;
  r[1] = dy;
  r[2] = dz;
  rn = sqrt(dx*dx+dy*dy+dz*dz);

  quatI = bonus[ellipsoid[i]].quat;
  quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];
  quatJ = bonus[ellipsoid[j]].quat;
  quatJtmp[0] = quatJ[0]; quatJtmp[1] = quatJ[1]; quatJtmp[2] = quatJ[2]; quatJtmp[3] = quatJ[3];

  LAMMPS_NS::quat_vec_rot(uI,u0,quatItmp);
  LAMMPS_NS::quat_vec_rot(uJ,u0,quatJtmp);

  phi = 0;
  theta1 = 0;
  theta2 = 0;

  // calc (and condition) a0, a1, a2 based on inputs
  if (rtype[m] == PHI) {
      a0 = LAMMPS_NS::vec_dot(uI,uJ);
      a1 = 0;
      a2 = 0;

      // calculate thetas and phi 
      if (a0 > 1) a0 = 1;
      if (a0 < -1) a0 = -1;
      if (a1 > 1) a1 = 1;
      if (a1 < -1) a1 = -1;
      if (a2 > 1) a2 = 1;
      if (a2 < -1) a2 = -1;
      phi = acos(a0);
      angles[m] = phi;
      dphi = phi - target[m];
      dtheta1 = 0;
      dtheta2 = 0;

  }
  else if (rtype[m] == THETA1) {
      a0 = 0;
      a1 = LAMMPS_NS::vec_dot(uI,r)/rn;
      a2 = 0;

      // note that thetas are pi minus the input value
      // this is since rij is defined as ri - rj (instead of the typical rj - ri)
      if (a0 > 1) a0 = 1;
      if (a0 < -1) a0 = -1;
      if (a1 > 1) a1 = 1;
      if (a1 < -1) a1 = -1;
      if (a2 > 1) a2 = 1;
      if (a2 < -1) a2 = -1;
      theta1 = acos(a1);
      angles[m] = theta1;
      dtheta1 = theta1 - target[m];
      dtheta2 = 0;
      dphi = 0;
  }
  else if (rtype[m] == THETA2) {
      a0 = 0;
      a1 = 0;
      a2 = LAMMPS_NS::vec_dot(uJ,r)/rn;
  
      if (a0 > 1) a0 = 1;
      if (a0 < -1) a0 = -1;
      if (a1 > 1) a1 = 1;
      if (a1 < -1) a1 = -1;
      if (a2 > 1) a2 = 1;
      if (a2 < -1) a2 = -1;
      theta2 = acos(a2);
      angles[m] = theta2;
      dtheta2 = theta2 - target[m];
      dtheta1 = 0;
      dphi = 0;
  }

  // ramp the spring constant 
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double kspr = kstart[m] + delta * (kstop[m] - kstart[m]);

  // calc energy
  U = 0.5*(kspr*dtheta1*dtheta1 + kspr*dtheta2*dtheta2 + kspr*dphi*dphi);

  // calc partial derivatives
  inverse = sqrt(1-a0*a0);
  if (inverse < SMALL) inverse = SMALL;
  dUda0 = - kspr*dphi / inverse;
  
  inverse = sqrt(1-a1*a1);
  if (inverse < SMALL) inverse = SMALL;
  dUda1 = - kspr*dtheta1 / inverse;
  
  inverse = sqrt(1-a2*a2);
  if (inverse < SMALL) inverse = SMALL;
  dUda2 = - kspr*dtheta2 / inverse;
  
  // forces 
  fterm1[0] = -(uI[0]/rn - r[0]*a1/rn/rn) * dUda1;
  fterm1[1] = -(uI[1]/rn - r[1]*a1/rn/rn) * dUda1;
  fterm1[2] = -(uI[2]/rn - r[2]*a1/rn/rn) * dUda1;
  
  fterm2[0] = -(uJ[0]/rn - r[0]*a2/rn/rn) * dUda2;
  fterm2[1] = -(uJ[1]/rn - r[1]*a2/rn/rn) * dUda2;
  fterm2[2] = -(uJ[2]/rn - r[2]*a2/rn/rn) * dUda2;
  
  force1[0] = fterm1[0] + fterm2[0];
  force1[1] = fterm1[1] + fterm2[1];
  force1[2] = fterm1[2] + fterm2[2];
  
  f[i][0] += force1[0];
  f[i][1] += force1[1];
  f[i][2] += force1[2];
  
  f[j][0] -= force1[0];
  f[j][1] -= force1[1];
  f[j][2] -= force1[2];
  
  //torque (see Allen and Tildsley Apendix C)
  dUdfI[0] = r[0] / rn * dUda1 + uJ[0]*dUda0;
  dUdfI[1] = r[1] / rn * dUda1 + uJ[1]*dUda0;
  dUdfI[2] = r[2] / rn * dUda1 + uJ[2]*dUda0;
  
  torqI[0] = -(uI[1]*dUdfI[2] - uI[2]*dUdfI[1]);
  torqI[1] = -(uI[2]*dUdfI[0] - uI[0]*dUdfI[2]);
  torqI[2] = -(uI[0]*dUdfI[1] - uI[1]*dUdfI[0]);
  
  dUdfJ[0] = r[0] / rn * dUda2 + uI[0]*dUda0;
  dUdfJ[1] = r[1] / rn * dUda2 + uI[1]*dUda0;
  dUdfJ[2] = r[2] / rn * dUda2 + uI[2]*dUda0;
  
  torqJ[0] = -(uJ[1]*dUdfJ[2] - uJ[2]*dUdfJ[1]);
  torqJ[1] = -(uJ[2]*dUdfJ[0] - uJ[0]*dUdfJ[2]);
  torqJ[2] = -(uJ[0]*dUdfJ[1] - uJ[1]*dUdfJ[0]);
  
  //when ii==0, apply to i; when ii==1 apply to j
  torque[i][0] += torqI[0];
  torque[i][1] += torqI[1];
  torque[i][2] += torqI[2];
  
  torque[j][0] += torqJ[0];
  torque[j][1] += torqJ[1];
  torque[j][2] += torqJ[2];
    
  energy+=U;
}


/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixRestrainOrient::compute_scalar()
{
  MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
  return energy_all;
}
/* ----------------------------------------------------------------------
   return angle of nth restraint
------------------------------------------------------------------------- */

double FixRestrainOrient::compute_vector(int n)
{
  return angles[n];
}
