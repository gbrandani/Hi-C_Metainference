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
   Contributing author: Joshua Lequieu (U Chicago)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_zewdie.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "integrate.h"
#include "citeme.h"
#include "memory.h"
#include "error.h"
#include "math_vector.h"

using namespace LAMMPS_NS;

static const char cite_pair_zewdie[] =
  "pair zewdie command:\n\n"
  "@Article{Zewdie1998,\n"
  " author =  {H. Zewdie},\n"
  " title =   {Computer-simulation studies of diskotic liquid crystals},\n"
  " journal = {Phys. Rev. E},\n"
  " year =    1998,\n"
  " volume =  57,\n"
  " pages =   {1793-1805}\n"
  "}\n\n";

/* ---------------------------------------------------------------------- */

PairZewdie::PairZewdie(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_zewdie);

  single_enable = 0;
  writedata = 1;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairZewdie::~PairZewdie()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(form);
    //memory->destroy(epsilon);
    //memory->destroy(sigma);
    memory->destroy(shape1);
    //memory->destroy(shape2);
    //memory->destroy(well);
    memory->destroy(cut);
    //memory->destroy(lj1);
    //memory->destroy(lj2);
    //memory->destroy(lj3);
    //memory->destroy(lj4);
    memory->destroy(offset);

    //delete zewdie parameters
    memory->destroy(epsilon0);
    memory->destroy(sigma0);

    //memory->destroy(sigma000);
    //memory->destroy(sigmacc2);
    //memory->destroy(sigma220);
    //memory->destroy(sigma222);
    //memory->destroy(sigma224);

    //memory->destroy(epsilon000);
    //memory->destroy(epsiloncc2);
    //memory->destroy(epsilon220);
    //memory->destroy(epsilon222);
    //memory->destroy(epsilon224);


    //delete [] lshape;
    //delete [] setwell;
  }
}

/* ---------------------------------------------------------------------- */

void PairZewdie::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double evdwl,one_eng,rsq,r2inv,r6inv,forcelj,factor_lj;
  double fforce[3],ttor[3],rtor[3],r12[3];
  //double a1[3][3],b1[3][3],g1[3][3],a2[3][3],b2[3][3],g2[3][3],temp[3][3];
  int *ilist,*jlist,*numneigh,**firstneigh;
  double *iquat,*jquat;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **f = atom->f;
  double **tor = atom->torque;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double sig,eps;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    if (form[itype][itype] == ELLIPSE_ELLIPSE) {
      iquat = bonus[ellipsoid[i]].quat;
      //MathExtra::quat_to_mat_trans(iquat,a1);
      //MathExtra::diag_times3(well[itype],a1,temp);
      //MathExtra::transpose_times3(a1,temp,b1);
      //MathExtra::diag_times3(shape2[itype],a1,temp);
      //MathExtra::transpose_times3(a1,temp,g1);
    }

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      // r12 = center to center vector

      r12[0] = x[j][0]-x[i][0];
      r12[1] = x[j][1]-x[i][1];
      r12[2] = x[j][2]-x[i][2];
      rsq = MathExtra::dot3(r12,r12);
      jtype = type[j];

      eps = epsilon0[itype][jtype];
      sig = sigma0[itype][jtype];
      
      // compute if less than cutoff

      if (rsq < cutsq[itype][jtype]) {

        switch (form[itype][jtype]) {
        case SPHERE_SPHERE:
          //r2inv = 1.0/rsq;
          //r6inv = r2inv*r2inv*r2inv;
          //forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          //forcelj *= -r2inv;
          //if (eflag) one_eng =
          //             r6inv*(r6inv*lj3[itype][jtype]-lj4[itype][jtype]) -
          //             offset[itype][jtype];
          //fforce[0] = r12[0]*forcelj;
          //fforce[1] = r12[1]*forcelj;
          //fforce[2] = r12[2]*forcelj;
          //ttor[0] = ttor[1] = ttor[2] = 0.0;
          //rtor[0] = rtor[1] = rtor[2] = 0.0;
          error->all(FLERR,"SPHERE_SPHERE interactions are not implemented");
          break;

        case SPHERE_ELLIPSE:
          jquat = bonus[ellipsoid[j]].quat;
          one_eng = zewdie_lj(i,j,jquat,r12,eps,sig,fforce,rtor);
          ttor[0] = ttor[1] = ttor[2] = 0.0;
          break;

        case ELLIPSE_SPHERE:
          one_eng = zewdie_lj(j,i,iquat,r12,eps,sig,fforce,ttor);
          rtor[0] = rtor[1] = rtor[2] = 0.0;
          break;

        default:
          jquat = bonus[ellipsoid[j]].quat;
          one_eng = zewdie_analytic(i,j,iquat,jquat,r12,eps,sig,fforce,ttor,rtor);
          break;
        }

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        ttor[0] *= factor_lj;
        ttor[1] *= factor_lj;
        ttor[2] *= factor_lj;

        f[i][0] += fforce[0];
        f[i][1] += fforce[1];
        f[i][2] += fforce[2];
        tor[i][0] += ttor[0];
        tor[i][1] += ttor[1];
        tor[i][2] += ttor[2];

        if (newton_pair || j < nlocal) {
          rtor[0] *= factor_lj;
          rtor[1] *= factor_lj;
          rtor[2] *= factor_lj;
          f[j][0] -= fforce[0];
          f[j][1] -= fforce[1];
          f[j][2] -= fforce[2];
          tor[j][0] += rtor[0];
          tor[j][1] += rtor[1];
          tor[j][2] += rtor[2];
        }

       if (eflag) evdwl = factor_lj*one_eng;

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                 -r12[0],-r12[1],-r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairZewdie::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(form,n+1,n+1,"pair:form");
  //memory->create(epsilon,n+1,n+1,"pair:epsilon");
  //memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(shape1,n+1,3,"pair:shape1");
  //memory->create(shape2,n+1,3,"pair:shape2");
  //memory->create(well,n+1,3,"pair:well");
  memory->create(cut,n+1,n+1,"pair:cut");
  //memory->create(lj1,n+1,n+1,"pair:lj1");
  //memory->create(lj2,n+1,n+1,"pair:lj2");
  //memory->create(lj3,n+1,n+1,"pair:lj3");
  //memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
  
  //zewdie parameters 
  // -set of \sigma0 and \epsilon0 parameters for EACH pairwise combination
  memory->create(epsilon0,n+1,n+1,"pair:epsilon0");
  memory->create(sigma0,n+1,n+1,"pair:sigma0");
  // - the other params will be const for the pair style
  //memory->create(epsilon000,n+1,n+1,"pair:epsilon000");
  //memory->create(epsiloncc2,n+1,n+1,"pair:epsiloncc2");
  //memory->create(epsilon220,n+1,n+1,"pair:epsilon220");
  //memory->create(epsilon222,n+1,n+1,"pair:epsilon222");
  //memory->create(epsilon224,n+1,n+1,"pair:epsilon224");
  //memory->create(sigma000,n+1,n+1,"pair:sigma000");
  //memory->create(sigmacc2,n+1,n+1,"pair:sigmacc2");
  //memory->create(sigma220,n+1,n+1,"pair:sigma220");
  //memory->create(sigma222,n+1,n+1,"pair:sigma222");
  //memory->create(sigma224,n+1,n+1,"pair:sigma224");

  //lshape = new double[n+1];
  //setwell = new int[n+1];
  //for (int i = 1; i <= n; i++) setwell[i] = 0;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairZewdie::settings(int narg, char **arg)
{
  if (narg != 10) error->all(FLERR,"Illegal pair_style command");

  //gamma = force->numeric(FLERR,arg[0]);
  //upsilon = force->numeric(FLERR,arg[1])/2.0;
  //mu = force->numeric(FLERR,arg[2]);
  epsilon000 = force->numeric(FLERR,arg[0]);
  epsiloncc2 = force->numeric(FLERR,arg[1]);
  epsilon220 = force->numeric(FLERR,arg[2]);
  epsilon222 = force->numeric(FLERR,arg[3]);
  epsilon224 = force->numeric(FLERR,arg[4]);

  sigma000   = force->numeric(FLERR,arg[5]);
  sigmacc2   = force->numeric(FLERR,arg[6]);
  sigma220   = force->numeric(FLERR,arg[7]);
  sigma222   = force->numeric(FLERR,arg[8]);
  sigma224   = force->numeric(FLERR,arg[9]);

  //cut_global = force->numeric(FLERR,arg[10]);

  // reset cutoffs that have been explicitly set

  //if (allocated) {
  //  int i,j;
  //  for (i = 1; i <= atom->ntypes; i++)
  //    for (j = i+1; j <= atom->ntypes; j++)
  //      if (setflag[i][j]) cut[i][j] = cut_global;
  //}
  //error->all(FLERR,"PairZewdie::function not checked");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairZewdie::coeff(int narg, char **arg)
{
  //if (narg < 4 || narg > 5)
  if (narg != 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  //double epsilon_one = force->numeric(FLERR,arg[2]);
  //double sigma_one = force->numeric(FLERR,arg[3]);
  //double eia_one = force->numeric(FLERR,arg[4]);
  //double eib_one = force->numeric(FLERR,arg[5]);
  //double eic_one = force->numeric(FLERR,arg[6]);
  //double eja_one = force->numeric(FLERR,arg[7]);
  //double ejb_one = force->numeric(FLERR,arg[8]);
  //double ejc_one = force->numeric(FLERR,arg[9]);

  double epsilon0_one   = force->numeric(FLERR,arg[2]);
  double sigma0_one     = force->numeric(FLERR,arg[3]);
  double cut_one        = force->numeric(FLERR,arg[4]);
  

  //if cut defined in coeffs, use it. Otherwise use the global cutoff
  //double cut_one = cut_global;
  //if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      //epsilon[i][j] = epsilon_one;
      //sigma[i][j] = sigma_one;
      //cut[i][j] = cut_one;
      epsilon0[i][j] = epsilon0_one;
      sigma0[i][j] = sigma0_one;
      cut[i][j] = cut_one;
  
      //if (eia_one != 0.0 || eib_one != 0.0 || eic_one != 0.0) {
      //  well[i][0] = pow(eia_one,-1.0/mu);
      //  well[i][1] = pow(eib_one,-1.0/mu);
      //  well[i][2] = pow(eic_one,-1.0/mu);
      //  if (eia_one == eib_one && eib_one == eic_one) setwell[i] = 2;
      //  else setwell[i] = 1;
      //}
      //if (eja_one != 0.0 || ejb_one != 0.0 || ejc_one != 0.0) {
      //  well[j][0] = pow(eja_one,-1.0/mu);
      //  well[j][1] = pow(ejb_one,-1.0/mu);
      //  well[j][2] = pow(ejc_one,-1.0/mu);
      //  if (eja_one == ejb_one && ejb_one == ejc_one) setwell[j] = 2;
      //  else setwell[j] = 1;
      //}

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
  
  //error->all(FLERR,"PairZewdie::function not checked");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairZewdie::init_style()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Pair zewdie requires atom style ellipsoid");

  neighbor->request(this,instance_me);

  // per-type shape precalculations
  // require that atom shapes are identical within each type
  // if shape = 0 for point particle, set shape = 1 as required by Gay-Berne

  for (int i = 1; i <= atom->ntypes; i++) {
    if (!atom->shape_consistency(i,shape1[i][0],shape1[i][1],shape1[i][2]))
      error->all(FLERR,
                 "Pair zewdie requires atoms with same type have same shape");

    //if its not a sphere, then check aspect ratio (expecting 1x2x2
    if ((shape1[i][0] != shape1[i][1]) || (shape1[i][0] != shape1[i][2])){
      if ((2*shape1[i][0] != shape1[i][1]) || 
          (2*shape1[i][0] != shape1[i][2]) || 
          (shape1[i][1] != shape1[i][2])){
          error->universe_warn(FLERR,"The zewdie usually has 2*xshape = yshape = zshape, but you don't. If you want to do this, you'll need different pair_style values for the new shape. The shape parameters dont effect the potential, so be careful!");
      }
    }
    //if (shape1[i][0] == 0.0)
    //  shape1[i][0] = shape1[i][1] = shape1[i][2] = 1.0;
    //shape2[i][0] = shape1[i][0]*shape1[i][0];
    //shape2[i][1] = shape1[i][1]*shape1[i][1];
    //shape2[i][2] = shape1[i][2]*shape1[i][2];
    //lshape[i] = (shape1[i][0]*shape1[i][1]+shape1[i][2]*shape1[i][2]) *
    //  sqrt(shape1[i][0]*shape1[i][1]);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairZewdie::init_one(int i, int j)
{
  //if (setwell[i] == 0 || setwell[j] == 0)
  //  error->all(FLERR,"Pair zewdie epsilon a,b,c coeffs are not all set");

  if (setflag[i][j] == 0) {
    epsilon0[i][j] = mix_energy(epsilon0[i][i],epsilon0[j][j],
                               sigma0[i][i],sigma0[j][j]);
    sigma0[i][j] = mix_distance(sigma0[i][i],sigma0[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  // perhaps reincorporate LJ precalculation again to speed compute evaluation
  // for now, I'll comment them out
  //lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  //lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  //lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  //lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  if (offset_flag) {
    //double ratio = sigma[i][j] / cut[i][j];
    //offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
    error->all(FLERR,"Offset flag set but not implemented.");
  } else offset[i][j] = 0.0;

  int ishape = 0;
  if (shape1[i][0] != shape1[i][1] ||
      shape1[i][0] != shape1[i][2] ||
      shape1[i][1] != shape1[i][2]) ishape = 1;
  //if (setwell[i] == 1) ishape = 1;
  int jshape = 0;
  if (shape1[j][0] != shape1[j][1] ||
      shape1[j][0] != shape1[j][2] ||
      shape1[j][1] != shape1[j][2]) jshape = 1;
  //if (setwell[j] == 1) jshape = 1;

  if (ishape == 0 && jshape == 0)
    form[i][i] = form[j][j] = form[i][j] = form[j][i] = SPHERE_SPHERE;
  else if (ishape == 0) {
    form[i][i] = SPHERE_SPHERE; form[j][j] = ELLIPSE_ELLIPSE;
    form[i][j] = SPHERE_ELLIPSE; form[j][i] = ELLIPSE_SPHERE;
  } else if (jshape == 0) {
    form[j][j] = SPHERE_SPHERE; form[i][i] = ELLIPSE_ELLIPSE;
    form[j][i] = SPHERE_ELLIPSE; form[i][j] = ELLIPSE_SPHERE;
  } else
    form[i][i] = form[j][j] = form[i][j] = form[j][i] = ELLIPSE_ELLIPSE;

  epsilon0[j][i] = epsilon0[i][j];
  sigma0[j][i] = sigma0[i][j];
  //lj1[j][i] = lj1[i][j];
  //lj2[j][i] = lj2[i][j];
  //lj3[j][i] = lj3[i][j];
  //lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairZewdie::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++) {
    //fwrite(&setwell[i],sizeof(int),1,fp);
    //if (setwell[i]) fwrite(&well[i][0],sizeof(double),3,fp);
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon0[i][j],sizeof(double),1,fp);
        fwrite(&sigma0[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairZewdie::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    //if (me == 0) fread(&setwell[i],sizeof(int),1,fp);
    //MPI_Bcast(&setwell[i],1,MPI_INT,0,world);
    //if (setwell[i]) {
    //  if (me == 0) fread(&well[i][0],sizeof(double),3,fp);
    //  MPI_Bcast(&well[i][0],3,MPI_DOUBLE,0,world);
    //}
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon0[i][j],sizeof(double),1,fp);
          fread(&sigma0[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
  }

}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairZewdie::write_restart_settings(FILE *fp)
{
  fwrite(&epsilon000,sizeof(double),1,fp);
  fwrite(&epsiloncc2,sizeof(double),1,fp);
  fwrite(&epsilon220,sizeof(double),1,fp);
  fwrite(&epsilon222,sizeof(double),1,fp);
  fwrite(&epsilon224,sizeof(double),1,fp);

  fwrite(&sigma000,sizeof(double),1,fp);
  fwrite(&sigmacc2,sizeof(double),1,fp);
  fwrite(&sigma220,sizeof(double),1,fp);
  fwrite(&sigma222,sizeof(double),1,fp);
  fwrite(&sigma224,sizeof(double),1,fp);
  
  //not sure if I need these 
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairZewdie::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&epsilon000,sizeof(double),1,fp);
    fread(&epsiloncc2,sizeof(double),1,fp);
    fread(&epsilon220,sizeof(double),1,fp);
    fread(&epsilon222,sizeof(double),1,fp);
    fread(&epsilon224,sizeof(double),1,fp);

    fread(&sigma000,sizeof(double),1,fp);
    fread(&sigmacc2,sizeof(double),1,fp);
    fread(&sigma220,sizeof(double),1,fp);
    fread(&sigma222,sizeof(double),1,fp);
    fread(&sigma224,sizeof(double),1,fp);
    
    //not sure if I need these 
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&epsilon000,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsiloncc2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilon220,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilon222,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&epsilon224,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&sigma000,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigmacc2,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma220,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma222,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&sigma224,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairZewdie::write_data(FILE *fp)
{
  error->universe_warn(FLERR,"zewdie write_data hasn't been tested, It SHOULD work though, but do some tests before trusting it.");

  for (int i = 1; i <= atom->ntypes; i++)
    //fprintf(fp,"%d %g %g %g %g %g %g %g %g\n",i,
    //        epsilon[i][i],sigma[i][i],
    //        pow(well[i][0],-mu),pow(well[i][1],-mu),pow(well[i][2],-mu),
    //        pow(well[i][0],-mu),pow(well[i][1],-mu),pow(well[i][2],-mu));
    fprintf(fp,"%d %g %g %g\n",i,
            epsilon0[i][i],sigma0[i][i],cut[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairZewdie::write_data_all(FILE *fp)
{
  error->universe_warn(FLERR,"zewdie write_data_all hasn't been tested, It SHOULD work though, but do some tests before trusting it.");
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      //fprintf(fp,"%d %d %g %g %g %g %g %g %g %g %g\n",i,j,
      //        epsilon[i][i],sigma[i][i],
      //        pow(well[i][0],-mu),pow(well[i][1],-mu),pow(well[i][2],-mu),
      //        pow(well[j][0],-mu),pow(well[j][1],-mu),pow(well[j][2],-mu),
      //        cut[i][j]);
      fprintf(fp,"%d %d %g %g %g\n",i,j,
              epsilon0[i][i],sigma0[i][i],cut[i][j]);
}

//josh's mods to the zewdie potential
//using forces and torques andres calculated
double PairZewdie::zewdie_analytic(const int i,const int j, double *quat1, 
				       double *quat2, double *r12,
                       const double epsilon, const double sigma, double *fforce,
                                       double *ttor, double *rtor)
{
  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;
  // Note: to accelerate this code, one could add statements like:
  //       "if (!(newton_pair || j < nlocal)) {}"
  //     to wrap code blocks that are unnecessary of newton_pair = false
  //     since some evaluations are not used.
  //     currently I dont do this, so newton_pair=false will give effectively no speedup .

  double potl;
  double uinit1[3],u1[3];
  double uinit2[3],u2[3];

  //assume that the orientation of the particles 
  uinit1[0] = 1; uinit1[1] = 0; uinit1[2] = 0; 
  uinit2[0] = 1; uinit2[1] = 0; uinit2[2] = 0;

  //move quat into a different data type so I can use quat_vec_rot
  double q1[4], q2[4];
  q1[0] = quat1[0]; q1[1] = quat1[1]; q1[2] = quat1[2]; q1[3] = quat1[3];
  q2[0] = quat2[0]; q2[1] = quat2[1]; q2[2] = quat2[2]; q2[3] = quat2[3];
  
  LAMMPS_NS::quat_vec_rot(u1,uinit1,q1);
  LAMMPS_NS::vec_norm(u1);
  LAMMPS_NS::quat_vec_rot(u2,uinit2,q2);
  LAMMPS_NS::vec_norm(u2);


  //Define parameters for Zewdie model
  double pe0, pe000, pecc2, pe220, pe222, pe224;
  double ps0, ps000, pscc2, ps220, ps222, ps224;

  ps0   = sigma;
  ps000 = sigma000;
  pscc2 = sigmacc2;
  ps220 = sigma220;
  ps222 = sigma222;
  ps224 = sigma224;
  pe0   = epsilon;
  pe000 = epsilon000;
  pecc2 = epsiloncc2;
  pe220 = epsilon220;
  pe222 = epsilon222;
  pe224 = epsilon224;
 
  double u1x, u1y, u1z, u2x, u2y, u2z; //local orientation vectors
  double f0,f1,f2;
  double sig,eps,U;
  double S000,S202,S022,S220,S222,S224; 
  double dS000df0,dS000df1,dS000df2,dS202df0,dS202df1,dS202df2,
         dS022df0,dS022df1,dS022df2,dS220df0,dS220df1,dS220df2,
         dS222df0,dS222df1,dS222df2,dS224df0,dS224df1,dS224df2;
  double dsdf0,dsdf1,dsdf2,dedf0,dedf1,dedf2,dUdrn,dUdf0,dUdf1,dUdf2;
  double fterm1[3],fterm2[3],fterm3[3],forc[3],torque[3];
  double r[3],rh[3],rn;
  double gradu1U[3],gradu2U[3];
 
  // Business time, calculate forces and torques

  u1x = u1[0]; u1y = u1[1]; u1z = u1[2];
  u2x = u2[0]; u2y = u2[1]; u2z = u2[2];
  r[0] = -r12[0]; 
  r[1] = -r12[1]; 
  r[2] = -r12[2];

  rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  rh[0] = r[0]/rn;
  rh[1] = r[1]/rn;
  rh[2] = r[2]/rn;
  f0 = u1x*u2x + u1y*u2y + u1z*u2z;
  f1 = u1x*rh[0] + u1y*rh[1] + u1z*rh[2];
  f2 = u2x*rh[0] + u2y*rh[1] + u2z*rh[2];
  if (f0 > 1) f0 = 1;
  if (f0 < -1) f0 = -1;
  if (f1 > 1) f1 = 1;
  if (f1 < -1) f1 = -1;
  if (f2 > 1) f2 = 1;
  if (f2 < -1) f2 = -1;

  S000 = 1.;
  S202 = (3*f1*f1 - 1)/(2.*sqrt(5));
  S022 = (3*f2*f2 - 1)/(2.*sqrt(5));
  S220 = (3*f0*f0 - 1)/(2.*sqrt(5));
  S222 = (2 - 3*f1*f1 - 3*f2*f2 - 3*f0*f0 + 9*f0*f1*f2)/sqrt(70);
  S224 = (1 + 2*f0*f0 - 5*f1*f1 - 5*f2*f2 - 20*f0*f1*f2 + 35*f1*f1*f2*f2)/(4.*sqrt(70));
  sig = ps0*(S000*ps000 + S220*ps220 + S222*ps222 + S224*ps224 + (S022 + S202)*pscc2);
  eps = pe0*(S000*pe000 + S220*pe220 + S222*pe222 + S224*pe224 + (S022 + S202)*pecc2);
  U = 4*eps*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));

  //derivatives, for forces and torques
  dS000df0 = 0;
  dS000df1 = 0;
  dS000df2 = 0;
  dS202df0 = 0;
  dS202df1 = (3*f1)/sqrt(5);
  dS202df2 = 0;
  dS022df0 = 0;
  dS022df1 = 0;
  dS022df2 = (3*f2)/sqrt(5);
  dS220df0 = (3*f0)/sqrt(5);
  dS220df1 = 0;
  dS220df2 = 0;
  dS222df0 = (-6*f0 + 9*f1*f2)/sqrt(70);
  dS222df1 = (-6*f1 + 9*f0*f2)/sqrt(70);
  dS222df2 = (9*f0*f1 - 6*f2)/sqrt(70);
  dS224df0 = (4*f0 - 20*f1*f2)/(4.*sqrt(70));
  dS224df1 = (-10*f1 - 20*f0*f2 + 70*f1*f2*f2)/(4.*sqrt(70));
  dS224df2 = (-20*f0*f1 - 10*f2 + 70*f1*f1*f2)/(4.*sqrt(70));

  dsdf0 = ps0*(dS000df0*ps000 + dS220df0*ps220 + dS222df0*ps222 
          + dS224df0*ps224 + (dS022df0 + dS202df0)*pscc2);
  dsdf1 = ps0*(dS000df1*ps000 + dS220df1*ps220 + dS222df1*ps222 
          + dS224df1*ps224 + (dS022df1 + dS202df1)*pscc2);
  dsdf2 = ps0*(dS000df2*ps000 + dS220df2*ps220 + dS222df2*ps222 
          + dS224df2*ps224 + (dS022df2 + dS202df2)*pscc2);
  dedf0 = pe0*(dS000df0*pe000 + dS220df0*pe220 + dS222df0*pe222 
          + dS224df0*pe224 + (dS022df0 + dS202df0)*pecc2);
  dedf1 = pe0*(dS000df1*pe000 + dS220df1*pe220 + dS222df1*pe222 
          + dS224df1*pe224 + (dS022df1 + dS202df1)*pecc2);
  dedf2 = pe0*(dS000df2*pe000 + dS220df2*pe220 + dS222df2*pe222 
          + dS224df2*pe224 + (dS022df2 + dS202df2)*pecc2);

  dUdrn = 4*eps*((-12*pow(ps0,12))/pow(ps0 + rn - sig,13) + (6*pow(ps0,6))/pow(ps0 + rn - sig,7));
  dUdf0 = 4*eps*((12*dsdf0*pow(ps0,12))/pow(ps0 + rn - sig,13) - (6*dsdf0*pow(ps0,6))/pow(ps0 + rn - sig,7)) 
          + 4*dedf0*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));
  dUdf1 = 4*eps*((12*dsdf1*pow(ps0,12))/pow(ps0 + rn - sig,13) - (6*dsdf1*pow(ps0,6))/pow(ps0 + rn - sig,7)) 
          + 4*dedf1*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));
  dUdf2 = 4*eps*((12*dsdf2*pow(ps0,12))/pow(ps0 + rn - sig,13) - (6*dsdf2*pow(ps0,6))/pow(ps0 + rn - sig,7)) 
          + 4*dedf2*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));

  //force (see Allen and Tildsley Apendix C)
  //f_{ij} = -f_{ji}, so only need to calculate force once
  fterm1[0] = -(r[0]/rn) * dUdrn;
  fterm1[1] = -(r[1]/rn) * dUdrn;
  fterm1[2] = -(r[2]/rn) * dUdrn;

  fterm2[0] = -(u1x/rn - r[0]*f1/rn/rn) * dUdf1;
  fterm2[1] = -(u1y/rn - r[1]*f1/rn/rn) * dUdf1;
  fterm2[2] = -(u1z/rn - r[2]*f1/rn/rn) * dUdf1;

  fterm3[0] = -(u2x/rn - r[0]*f2/rn/rn) * dUdf2;
  fterm3[1] = -(u2y/rn - r[1]*f2/rn/rn) * dUdf2;
  fterm3[2] = -(u2z/rn - r[2]*f2/rn/rn) * dUdf2;

  forc[0] = fterm1[0] + fterm2[0] + fterm3[0];
  forc[1] = fterm1[1] + fterm2[1] + fterm3[1];
  forc[2] = fterm1[2] + fterm2[2] + fterm3[2];
  
  fforce[0] = forc[0]; 
  fforce[1] = forc[1];
  fforce[2] = forc[2]; 

  //torque (see Allen and Tildsley Apendix C)
  // tau_{ij} != -tau_{ji}, so need to calculate twice
  gradu1U[0] = (r[0]/ rn * dUdf1 + u2x*dUdf0);
  gradu1U[1] = (r[1]/ rn * dUdf1 + u2y*dUdf0);
  gradu1U[2] = (r[2]/ rn * dUdf1 + u2z*dUdf0);

  ttor[0] = -(u1y*gradu1U[2] - u1z*gradu1U[1]);
  ttor[1] = -(u1z*gradu1U[0] - u1x*gradu1U[2]);
  ttor[2] = -(u1x*gradu1U[1] - u1y*gradu1U[0]);

  gradu2U[0] = (r[0]/ rn * dUdf2 + u1x*dUdf0);
  gradu2U[1] = (r[1]/ rn * dUdf2 + u1y*dUdf0);
  gradu2U[2] = (r[2]/ rn * dUdf2 + u1z*dUdf0);

  rtor[0] = -(u2y*gradu2U[2] - u2z*gradu2U[1]);
  rtor[1] = -(u2z*gradu2U[0] - u2x*gradu2U[2]);
  rtor[2] = -(u2x*gradu2U[1] - u2y*gradu2U[0]);

  //return energy
  return U;
}


/* ----------------------------------------------------------------------
   compute analytic energy, force (fforce), and torque (rtor)
   between zewdie particle and lj particle

   here I assume that "i" is the lj particle and "j" is the zewdie particle
------------------------------------------------------------------------- */
double PairZewdie::zewdie_lj(const int i,const int j, 
				       double *quat2, double *r12,
                       const double epsilon, const double sigma, double *fforce,
                       double *rtor)
{
  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;

  // Note: to accelerate this code, one could add statements like:
  //       "if (!(newton_pair || j < nlocal)) {}"
  //     to wrap code blocks that are unnecessary of newton_pair = false
  //     since some evaluations are not used.
  //     currently I dont do this, so newton_pair=false will give effectively no speedup .

  double potl;
  double uinit2[3],u2[3];

  //assume that the orientation of the particles 
  uinit2[0] = 1; uinit2[1] = 0; uinit2[2] = 0;

  //move quat into a different data type so I can use quat_vec_rot
  double q2[4];
  q2[0] = quat2[0]; q2[1] = quat2[1]; q2[2] = quat2[2]; q2[3] = quat2[3];
  
  LAMMPS_NS::quat_vec_rot(u2,uinit2,q2);
  LAMMPS_NS::vec_norm(u2);


  //Define parameters for Zewdie model
  double pe0, pe000, pecc2, pe220, pe222, pe224;
  double ps0, ps000, pscc2, ps220, ps222, ps224;

  ps0   = sigma;
  ps000 = sigma000;
  pscc2 = sigmacc2;
  ps220 = sigma220;
  ps222 = sigma222;
  ps224 = sigma224;
  pe0   = epsilon;
  pe000 = epsilon000;
  pecc2 = epsiloncc2;
  pe220 = epsilon220;
  pe222 = epsilon222;
  pe224 = epsilon224;
 
  double u1x, u1y, u1z, u2x, u2y, u2z; //local orientation vectors
  double f0,f1,f2;
  double sig,eps,U;
  double S000,S202,S022,S220,S222,S224; 
  double dS000df0,dS000df1,dS000df2,dS202df0,dS202df1,dS202df2,
         dS022df0,dS022df1,dS022df2,dS220df0,dS220df1,dS220df2,
         dS222df0,dS222df1,dS222df2,dS224df0,dS224df1,dS224df2;
  double dsdf0,dsdf1,dsdf2,dedf0,dedf1,dedf2,dUdrn,dUdf0,dUdf1,dUdf2;
  double fterm1[3],fterm2[3],fterm3[3],forc[3],torque[3];
  double r[3],rh[3],rn;
  double gradu1U[3],gradu2U[3];
 
  // Business time, calculate forces and torques

  // The big change for the sphere_zewdie interaction:
  // u1 == rh, a0 == a2

  u2x = u2[0]; u2y = u2[1]; u2z = u2[2];
  r[0] = -r12[0]; 
  r[1] = -r12[1]; 
  r[2] = -r12[2];

  rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  rh[0] = r[0]/rn;
  rh[1] = r[1]/rn;
  rh[2] = r[2]/rn;

  u1x = rh[0]; u1y = rh[1]; u1z = rh[2];

  //f0 = u1x*u2x + u1y*u2y + u1z*u2z;
  f1 = u1x*rh[0] + u1y*rh[1] + u1z*rh[2];
  f2 = u2x*rh[0] + u2y*rh[1] + u2z*rh[2];
  f0 = f2;
  if (f0 > 1) f0 = 1;
  if (f0 < -1) f0 = -1;
  if (f1 > 1) f1 = 1;
  if (f1 < -1) f1 = -1;
  if (f2 > 1) f2 = 1;
  if (f2 < -1) f2 = -1;

  S000 = 1.;
  S202 = (3*f1*f1 - 1)/(2.*sqrt(5));
  S022 = (3*f2*f2 - 1)/(2.*sqrt(5));
  S220 = (3*f0*f0 - 1)/(2.*sqrt(5));
  S222 = (2 - 3*f1*f1 - 3*f2*f2 - 3*f0*f0 + 9*f0*f1*f2)/sqrt(70);
  S224 = (1 + 2*f0*f0 - 5*f1*f1 - 5*f2*f2 - 20*f0*f1*f2 + 35*f1*f1*f2*f2)/(4.*sqrt(70));
  sig = ps0*(S000*ps000 + S220*ps220 + S222*ps222 + S224*ps224 + (S022 + S202)*pscc2);
  eps = pe0*(S000*pe000 + S220*pe220 + S222*pe222 + S224*pe224 + (S022 + S202)*pecc2);
  U = 4*eps*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));

  //derivatives, for forces and torques
  dS000df0 = 0;
  dS000df1 = 0;
  dS000df2 = 0;
  dS202df0 = 0;
  dS202df1 = (3*f1)/sqrt(5);
  dS202df2 = 0;
  dS022df0 = 0;
  dS022df1 = 0;
  dS022df2 = (3*f2)/sqrt(5);
  dS220df0 = (3*f0)/sqrt(5);
  dS220df1 = 0;
  dS220df2 = 0;
  dS222df0 = (-6*f0 + 9*f1*f2)/sqrt(70);
  dS222df1 = (-6*f1 + 9*f0*f2)/sqrt(70);
  dS222df2 = (9*f0*f1 - 6*f2)/sqrt(70);
  dS224df0 = (4*f0 - 20*f1*f2)/(4.*sqrt(70));
  dS224df1 = (-10*f1 - 20*f0*f2 + 70*f1*f2*f2)/(4.*sqrt(70));
  dS224df2 = (-20*f0*f1 - 10*f2 + 70*f1*f1*f2)/(4.*sqrt(70));

  dsdf0 = ps0*(dS000df0*ps000 + dS220df0*ps220 + dS222df0*ps222 
          + dS224df0*ps224 + (dS022df0 + dS202df0)*pscc2);
  dsdf1 = ps0*(dS000df1*ps000 + dS220df1*ps220 + dS222df1*ps222 
          + dS224df1*ps224 + (dS022df1 + dS202df1)*pscc2);
  dsdf2 = ps0*(dS000df2*ps000 + dS220df2*ps220 + dS222df2*ps222 
          + dS224df2*ps224 + (dS022df2 + dS202df2)*pscc2);
  dedf0 = pe0*(dS000df0*pe000 + dS220df0*pe220 + dS222df0*pe222 
          + dS224df0*pe224 + (dS022df0 + dS202df0)*pecc2);
  dedf1 = pe0*(dS000df1*pe000 + dS220df1*pe220 + dS222df1*pe222 
          + dS224df1*pe224 + (dS022df1 + dS202df1)*pecc2);
  dedf2 = pe0*(dS000df2*pe000 + dS220df2*pe220 + dS222df2*pe222 
          + dS224df2*pe224 + (dS022df2 + dS202df2)*pecc2);

  dUdrn = 4*eps*((-12*pow(ps0,12))/pow(ps0 + rn - sig,13) + (6*pow(ps0,6))/pow(ps0 + rn - sig,7));
  dUdf0 = 4*eps*((12*dsdf0*pow(ps0,12))/pow(ps0 + rn - sig,13) - (6*dsdf0*pow(ps0,6))/pow(ps0 + rn - sig,7)) 
          + 4*dedf0*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));
  dUdf1 = 4*eps*((12*dsdf1*pow(ps0,12))/pow(ps0 + rn - sig,13) - (6*dsdf1*pow(ps0,6))/pow(ps0 + rn - sig,7)) 
          + 4*dedf1*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));
  dUdf2 = 4*eps*((12*dsdf2*pow(ps0,12))/pow(ps0 + rn - sig,13) - (6*dsdf2*pow(ps0,6))/pow(ps0 + rn - sig,7)) 
          + 4*dedf2*(pow(ps0,12)/pow(ps0 + rn - sig,12) - pow(ps0,6)/pow(ps0 + rn - sig,6));

  //force (see Allen and Tildsley Apendix C)
  //f_{ij} = -f_{ji}, so only need to calculate force once
  fterm1[0] = -(r[0]/rn) * dUdrn;
  fterm1[1] = -(r[1]/rn) * dUdrn;
  fterm1[2] = -(r[2]/rn) * dUdrn;

  fterm2[0] = -(u1x/rn - r[0]*f1/rn/rn) * dUdf1;
  fterm2[1] = -(u1y/rn - r[1]*f1/rn/rn) * dUdf1;
  fterm2[2] = -(u1z/rn - r[2]*f1/rn/rn) * dUdf1;

  fterm3[0] = -(u2x/rn - r[0]*f2/rn/rn) * dUdf2;
  fterm3[1] = -(u2y/rn - r[1]*f2/rn/rn) * dUdf2;
  fterm3[2] = -(u2z/rn - r[2]*f2/rn/rn) * dUdf2;

  forc[0] = fterm1[0] + fterm2[0] + fterm3[0];
  forc[1] = fterm1[1] + fterm2[1] + fterm3[1];
  forc[2] = fterm1[2] + fterm2[2] + fterm3[2];
  
  fforce[0] = forc[0]; 
  fforce[1] = forc[1];
  fforce[2] = forc[2]; 

  //torque (see Allen and Tildsley Apendix C)
  // tau_{ij} != -tau_{ji}, so need to calculate twice
  gradu1U[0] = (r[0]/ rn * dUdf1 + u2x*dUdf0);
  gradu1U[1] = (r[1]/ rn * dUdf1 + u2y*dUdf0);
  gradu1U[2] = (r[2]/ rn * dUdf1 + u2z*dUdf0);

  //ttor should be zero, but calculate it just in case (for debugging)
  double ttor[3];
  ttor[0] = -(u1y*gradu1U[2] - u1z*gradu1U[1]);
  ttor[1] = -(u1z*gradu1U[0] - u1x*gradu1U[2]);
  ttor[2] = -(u1x*gradu1U[1] - u1y*gradu1U[0]);

  gradu2U[0] = (r[0]/ rn * dUdf2 + u1x*dUdf0);
  gradu2U[1] = (r[1]/ rn * dUdf2 + u1y*dUdf0);
  gradu2U[2] = (r[2]/ rn * dUdf2 + u1z*dUdf0);

  rtor[0] = -(u2y*gradu2U[2] - u2z*gradu2U[1]);
  rtor[1] = -(u2z*gradu2U[0] - u2x*gradu2U[2]);
  rtor[2] = -(u2x*gradu2U[1] - u2y*gradu2U[0]);

  //return energy
  return U;
}

////original gayberne, not used
///* ----------------------------------------------------------------------
//   compute analytic energy, force (fforce), and torque (ttor & rtor)
//   based on rotation matrices a and precomputed matrices b and g
//   if newton is off, rtor is not calculated for ghost atoms
//------------------------------------------------------------------------- */
//
//double PairZewdie::zewdie_analytic(const int i,const int j,double a1[3][3],
//                                       double a2[3][3], double b1[3][3],
//                                       double b2[3][3], double g1[3][3],
//                                       double g2[3][3], double *r12,
//                                       const double rsq, double *fforce,
//                                       double *ttor, double *rtor)
//{
//  double tempv[3], tempv2[3];
//  double temp[3][3];
//  double temp1,temp2,temp3;
//
//  int *type = atom->type;
//  int newton_pair = force->newton_pair;
//  int nlocal = atom->nlocal;
//
//  double r12hat[3];
//  MathExtra::normalize3(r12,r12hat);
//  double r = sqrt(rsq);
//
//  // compute distance of closest approach
//
//  double g12[3][3];
//  MathExtra::plus3(g1,g2,g12);
//  double kappa[3];
//  int ierror = MathExtra::mldivide3(g12,r12,kappa);
//  if (ierror) error->all(FLERR,"Bad matrix inversion in mldivide3");
//
//  // tempv = G12^-1*r12hat
//
//  tempv[0] = kappa[0]/r;
//  tempv[1] = kappa[1]/r;
//  tempv[2] = kappa[2]/r;
//  double sigma12 = MathExtra::dot3(r12hat,tempv);
//  sigma12 = pow(0.5*sigma12,-0.5);
//  double h12 = r-sigma12;
//
//  // energy
//  // compute u_r
//
//  double varrho = sigma[type[i]][type[j]]/(h12+gamma*sigma[type[i]][type[j]]);
//  double varrho6 = pow(varrho,6.0);
//  double varrho12 = varrho6*varrho6;
//  double u_r = 4.0*epsilon[type[i]][type[j]]*(varrho12-varrho6);
//
//  // compute eta_12
//
//  double eta = 2.0*lshape[type[i]]*lshape[type[j]];
//  double det_g12 = MathExtra::det3(g12);
//  eta = pow(eta/det_g12,upsilon);
//
//  // compute chi_12
//
//  double b12[3][3];
//  double iota[3];
//  MathExtra::plus3(b1,b2,b12);
//  ierror = MathExtra::mldivide3(b12,r12,iota);
//  if (ierror) error->all(FLERR,"Bad matrix inversion in mldivide3");
//
//  // tempv = G12^-1*r12hat
//
//  tempv[0] = iota[0]/r;
//  tempv[1] = iota[1]/r;
//  tempv[2] = iota[2]/r;
//  double chi = MathExtra::dot3(r12hat,tempv);
//  chi = pow(chi*2.0,mu);
//
//  // force
//  // compute dUr/dr
//
//  temp1 = (2.0*varrho12*varrho-varrho6*varrho)/sigma[type[i]][type[j]];
//  temp1 = temp1*24.0*epsilon[type[i]][type[j]];
//  double u_slj = temp1*pow(sigma12,3.0)/2.0;
//  double dUr[3];
//  temp2 = MathExtra::dot3(kappa,r12hat);
//  double uslj_rsq = u_slj/rsq;
//  dUr[0] = temp1*r12hat[0]+uslj_rsq*(kappa[0]-temp2*r12hat[0]);
//  dUr[1] = temp1*r12hat[1]+uslj_rsq*(kappa[1]-temp2*r12hat[1]);
//  dUr[2] = temp1*r12hat[2]+uslj_rsq*(kappa[2]-temp2*r12hat[2]);
//
//  // compute dChi_12/dr
//
//  double dchi[3];
//  temp1 = MathExtra::dot3(iota,r12hat);
//  temp2 = -4.0/rsq*mu*pow(chi,(mu-1.0)/mu);
//  dchi[0] = temp2*(iota[0]-temp1*r12hat[0]);
//  dchi[1] = temp2*(iota[1]-temp1*r12hat[1]);
//  dchi[2] = temp2*(iota[2]-temp1*r12hat[2]);
//
//  temp1 = -eta*u_r;
//  temp2 = eta*chi;
//  fforce[0] = temp1*dchi[0]-temp2*dUr[0];
//  fforce[1] = temp1*dchi[1]-temp2*dUr[1];
//  fforce[2] = temp1*dchi[2]-temp2*dUr[2];
//
//  // torque for particle 1 and 2
//  // compute dUr
//
//  tempv[0] = -uslj_rsq*kappa[0];
//  tempv[1] = -uslj_rsq*kappa[1];
//  tempv[2] = -uslj_rsq*kappa[2];
//  MathExtra::vecmat(kappa,g1,tempv2);
//  MathExtra::cross3(tempv,tempv2,dUr);
//  double dUr2[3];
//
//  if (newton_pair || j < nlocal) {
//    MathExtra::vecmat(kappa,g2,tempv2);
//    MathExtra::cross3(tempv,tempv2,dUr2);
//  }
//
//  // compute d_chi
//
//  MathExtra::vecmat(iota,b1,tempv);
//  MathExtra::cross3(tempv,iota,dchi);
//  temp1 = -4.0/rsq;
//  dchi[0] *= temp1;
//  dchi[1] *= temp1;
//  dchi[2] *= temp1;
//  double dchi2[3];
//
//  if (newton_pair || j < nlocal) {
//    MathExtra::vecmat(iota,b2,tempv);
//    MathExtra::cross3(tempv,iota,dchi2);
//    dchi2[0] *= temp1;
//    dchi2[1] *= temp1;
//    dchi2[2] *= temp1;
//  }
//
//  // compute d_eta
//
//  double deta[3];
//  deta[0] = deta[1] = deta[2] = 0.0;
//  compute_eta_torque(g12,a1,shape2[type[i]],temp);
//  temp1 = -eta*upsilon;
//  for (int m = 0; m < 3; m++) {
//    for (int y = 0; y < 3; y++) tempv[y] = temp1*temp[m][y];
//    MathExtra::cross3(a1[m],tempv,tempv2);
//    deta[0] += tempv2[0];
//    deta[1] += tempv2[1];
//    deta[2] += tempv2[2];
//  }
//
//  // compute d_eta for particle 2
//
//  double deta2[3];
//  if (newton_pair || j < nlocal) {
//    deta2[0] = deta2[1] = deta2[2] = 0.0;
//    compute_eta_torque(g12,a2,shape2[type[j]],temp);
//    for (int m = 0; m < 3; m++) {
//      for (int y = 0; y < 3; y++) tempv[y] = temp1*temp[m][y];
//      MathExtra::cross3(a2[m],tempv,tempv2);
//      deta2[0] += tempv2[0];
//      deta2[1] += tempv2[1];
//      deta2[2] += tempv2[2];
//    }
//  }
//
//  // torque
//
//  temp1 = u_r*eta;
//  temp2 = u_r*chi;
//  temp3 = chi*eta;
//
//  ttor[0] = (temp1*dchi[0]+temp2*deta[0]+temp3*dUr[0]) * -1.0;
//  ttor[1] = (temp1*dchi[1]+temp2*deta[1]+temp3*dUr[1]) * -1.0;
//  ttor[2] = (temp1*dchi[2]+temp2*deta[2]+temp3*dUr[2]) * -1.0;
//
//  if (newton_pair || j < nlocal) {
//    rtor[0] = (temp1*dchi2[0]+temp2*deta2[0]+temp3*dUr2[0]) * -1.0;
//    rtor[1] = (temp1*dchi2[1]+temp2*deta2[1]+temp3*dUr2[1]) * -1.0;
//    rtor[2] = (temp1*dchi2[2]+temp2*deta2[2]+temp3*dUr2[2]) * -1.0;
//  }
//
//  return temp1*chi;
//}


///* ----------------------------------------------------------------------
//   compute analytic energy, force (fforce), and torque (ttor)
//   between ellipsoid and lj particle
//------------------------------------------------------------------------- */
//
//double PairZewdie::zewdie_lj(const int i,const int j,double a1[3][3],
//                                 double b1[3][3],double g1[3][3],
//                                 double *r12,const double rsq,double *fforce,
//                                 double *ttor)
//{
//  //error->all(FLERR,"PairZewdie::zeqdie_lj not implemented");
//
//  double tempv[3], tempv2[3];
//  double temp[3][3];
//  double temp1,temp2,temp3;
//
//  int *type = atom->type;
//
//  double r12hat[3];
//  MathExtra::normalize3(r12,r12hat);
//  double r = sqrt(rsq);
//
//  // compute distance of closest approach
//
//  double g12[3][3];
//  g12[0][0] = g1[0][0]+shape2[type[j]][0];
//  g12[1][1] = g1[1][1]+shape2[type[j]][0];
//  g12[2][2] = g1[2][2]+shape2[type[j]][0];
//  g12[0][1] = g1[0][1]; g12[1][0] = g1[1][0];
//  g12[0][2] = g1[0][2]; g12[2][0] = g1[2][0];
//  g12[1][2] = g1[1][2]; g12[2][1] = g1[2][1];
//  double kappa[3];
//  int ierror = MathExtra::mldivide3(g12,r12,kappa);
//  if (ierror) error->all(FLERR,"Bad matrix inversion in mldivide3");
//
//  // tempv = G12^-1*r12hat
//
//  tempv[0] = kappa[0]/r;
//  tempv[1] = kappa[1]/r;
//  tempv[2] = kappa[2]/r;
//  double sigma12 = MathExtra::dot3(r12hat,tempv);
//  sigma12 = pow(0.5*sigma12,-0.5);
//  double h12 = r-sigma12;
//
//  // energy
//  // compute u_r
//
//  double varrho = sigma[type[i]][type[j]]/(h12+gamma*sigma[type[i]][type[j]]);
//  double varrho6 = pow(varrho,6.0);
//  double varrho12 = varrho6*varrho6;
//  double u_r = 4.0*epsilon[type[i]][type[j]]*(varrho12-varrho6);
//
//  // compute eta_12
//
//  double eta = 2.0*lshape[type[i]]*lshape[type[j]];
//  double det_g12 = MathExtra::det3(g12);
//  eta = pow(eta/det_g12,upsilon);
//
//  // compute chi_12
//
//  double b12[3][3];
//  double iota[3];
//  b12[0][0] = b1[0][0] + well[type[j]][0];
//  b12[1][1] = b1[1][1] + well[type[j]][0];
//  b12[2][2] = b1[2][2] + well[type[j]][0];
//  b12[0][1] = b1[0][1]; b12[1][0] = b1[1][0];
//  b12[0][2] = b1[0][2]; b12[2][0] = b1[2][0];
//  b12[1][2] = b1[1][2]; b12[2][1] = b1[2][1];
//  ierror = MathExtra::mldivide3(b12,r12,iota);
//  if (ierror) error->all(FLERR,"Bad matrix inversion in mldivide3");
//
//  // tempv = G12^-1*r12hat
//
//  tempv[0] = iota[0]/r;
//  tempv[1] = iota[1]/r;
//  tempv[2] = iota[2]/r;
//  double chi = MathExtra::dot3(r12hat,tempv);
//  chi = pow(chi*2.0,mu);
//
//  // force
//  // compute dUr/dr
//
//  temp1 = (2.0*varrho12*varrho-varrho6*varrho)/sigma[type[i]][type[j]];
//  temp1 = temp1*24.0*epsilon[type[i]][type[j]];
//  double u_slj = temp1*pow(sigma12,3.0)/2.0;
//  double dUr[3];
//  temp2 = MathExtra::dot3(kappa,r12hat);
//  double uslj_rsq = u_slj/rsq;
//  dUr[0] = temp1*r12hat[0]+uslj_rsq*(kappa[0]-temp2*r12hat[0]);
//  dUr[1] = temp1*r12hat[1]+uslj_rsq*(kappa[1]-temp2*r12hat[1]);
//  dUr[2] = temp1*r12hat[2]+uslj_rsq*(kappa[2]-temp2*r12hat[2]);
//
//  // compute dChi_12/dr
//
//  double dchi[3];
//  temp1 = MathExtra::dot3(iota,r12hat);
//  temp2 = -4.0/rsq*mu*pow(chi,(mu-1.0)/mu);
//  dchi[0] = temp2*(iota[0]-temp1*r12hat[0]);
//  dchi[1] = temp2*(iota[1]-temp1*r12hat[1]);
//  dchi[2] = temp2*(iota[2]-temp1*r12hat[2]);
//
//  temp1 = -eta*u_r;
//  temp2 = eta*chi;
//  fforce[0] = temp1*dchi[0]-temp2*dUr[0];
//  fforce[1] = temp1*dchi[1]-temp2*dUr[1];
//  fforce[2] = temp1*dchi[2]-temp2*dUr[2];
//
//  // torque for particle 1 and 2
//  // compute dUr
//
//  tempv[0] = -uslj_rsq*kappa[0];
//  tempv[1] = -uslj_rsq*kappa[1];
//  tempv[2] = -uslj_rsq*kappa[2];
//  MathExtra::vecmat(kappa,g1,tempv2);
//  MathExtra::cross3(tempv,tempv2,dUr);
//
//  // compute d_chi
//
//  MathExtra::vecmat(iota,b1,tempv);
//  MathExtra::cross3(tempv,iota,dchi);
//  temp1 = -4.0/rsq;
//  dchi[0] *= temp1;
//  dchi[1] *= temp1;
//  dchi[2] *= temp1;
//
//  // compute d_eta
//
//  double deta[3];
//  deta[0] = deta[1] = deta[2] = 0.0;
//  compute_eta_torque(g12,a1,shape2[type[i]],temp);
//  temp1 = -eta*upsilon;
//  for (int m = 0; m < 3; m++) {
//    for (int y = 0; y < 3; y++) tempv[y] = temp1*temp[m][y];
//    MathExtra::cross3(a1[m],tempv,tempv2);
//    deta[0] += tempv2[0];
//    deta[1] += tempv2[1];
//    deta[2] += tempv2[2];
//  }
//
//  // torque
//
//  temp1 = u_r*eta;
//  temp2 = u_r*chi;
//  temp3 = chi*eta;
//
//  ttor[0] = (temp1*dchi[0]+temp2*deta[0]+temp3*dUr[0]) * -1.0;
//  ttor[1] = (temp1*dchi[1]+temp2*deta[1]+temp3*dUr[1]) * -1.0;
//  ttor[2] = (temp1*dchi[2]+temp2*deta[2]+temp3*dUr[2]) * -1.0;
//
//  return temp1*chi;
//}
//
//
//// From gayberne, not needed for zewdie
///* ----------------------------------------------------------------------
//   torque contribution from eta
//   computes trace in the last doc equation for the torque derivative
//   code comes from symbolic solver dump
//   m is g12, m2 is a_i, s is the shape for the particle
//------------------------------------------------------------------------- */
//
//void PairZewdie::compute_eta_torque(double m[3][3], double m2[3][3],
//                                      double *s, double ans[3][3])
//{
//  double den = m[1][0]*m[0][2]*m[2][1]-m[0][0]*m[1][2]*m[2][1]-
//    m[0][2]*m[2][0]*m[1][1]+m[0][1]*m[2][0]*m[1][2]-
//    m[1][0]*m[0][1]*m[2][2]+m[0][0]*m[1][1]*m[2][2];
//
//  ans[0][0] = s[0]*(m[1][2]*m[0][1]*m2[0][2]+2.0*m[1][1]*m[2][2]*m2[0][0]-
//                    m[1][1]*m2[0][2]*m[0][2]-2.0*m[1][2]*m2[0][0]*m[2][1]+
//                    m2[0][1]*m[0][2]*m[2][1]-m2[0][1]*m[0][1]*m[2][2]-
//                    m[1][0]*m[2][2]*m2[0][1]+m[2][0]*m[1][2]*m2[0][1]+
//                    m[1][0]*m2[0][2]*m[2][1]-m2[0][2]*m[2][0]*m[1][1])/den;
//
//  ans[0][1] = s[0]*(m[0][2]*m2[0][0]*m[2][1]-m[2][2]*m2[0][0]*m[0][1]+
//                    2.0*m[0][0]*m[2][2]*m2[0][1]-m[0][0]*m2[0][2]*m[1][2]-
//                    2.0*m[2][0]*m[0][2]*m2[0][1]+m2[0][2]*m[1][0]*m[0][2]-
//                    m[2][2]*m[1][0]*m2[0][0]+m[2][0]*m2[0][0]*m[1][2]+
//                    m[2][0]*m2[0][2]*m[0][1]-m2[0][2]*m[0][0]*m[2][1])/den;
//
//  ans[0][2] = s[0]*(m[0][1]*m[1][2]*m2[0][0]-m[0][2]*m2[0][0]*m[1][1]-
//                    m[0][0]*m[1][2]*m2[0][1]+m[1][0]*m[0][2]*m2[0][1]-
//                    m2[0][1]*m[0][0]*m[2][1]-m[2][0]*m[1][1]*m2[0][0]+
//                    2.0*m[1][1]*m[0][0]*m2[0][2]-2.0*m[1][0]*m2[0][2]*m[0][1]+
//                    m[1][0]*m[2][1]*m2[0][0]+m[2][0]*m2[0][1]*m[0][1])/den;
//
//  ans[1][0] = s[1]*(-m[1][1]*m2[1][2]*m[0][2]+2.0*m[1][1]*m[2][2]*m2[1][0]+
//                    m[1][2]*m[0][1]*m2[1][2]-2.0*m[1][2]*m2[1][0]*m[2][1]+
//                    m2[1][1]*m[0][2]*m[2][1]-m2[1][1]*m[0][1]*m[2][2]-
//                    m[1][0]*m[2][2]*m2[1][1]+m[2][0]*m[1][2]*m2[1][1]-
//                    m2[1][2]*m[2][0]*m[1][1]+m[1][0]*m2[1][2]*m[2][1])/den;
//
//  ans[1][1] = s[1]*(m[0][2]*m2[1][0]*m[2][1]-m[0][1]*m[2][2]*m2[1][0]+
//                    2.0*m[2][2]*m[0][0]*m2[1][1]-m2[1][2]*m[0][0]*m[1][2]-
//                    2.0*m[2][0]*m2[1][1]*m[0][2]-m[1][0]*m[2][2]*m2[1][0]+
//                    m[2][0]*m[1][2]*m2[1][0]+m[1][0]*m2[1][2]*m[0][2]-
//                    m[0][0]*m2[1][2]*m[2][1]+m2[1][2]*m[0][1]*m[2][0])/den;
//
//  ans[1][2] = s[1]*(m[0][1]*m[1][2]*m2[1][0]-m[0][2]*m2[1][0]*m[1][1]-
//                    m[0][0]*m[1][2]*m2[1][1]+m[1][0]*m[0][2]*m2[1][1]+
//                    2.0*m[1][1]*m[0][0]*m2[1][2]-m[0][0]*m2[1][1]*m[2][1]+
//                    m[0][1]*m[2][0]*m2[1][1]-m2[1][0]*m[2][0]*m[1][1]-
//                    2.0*m[1][0]*m[0][1]*m2[1][2]+m[1][0]*m2[1][0]*m[2][1])/den;
//
//  ans[2][0] = s[2]*(-m[1][1]*m[0][2]*m2[2][2]+m[0][1]*m[1][2]*m2[2][2]+
//                    2.0*m[1][1]*m2[2][0]*m[2][2]-m[0][1]*m2[2][1]*m[2][2]+
//                    m[0][2]*m[2][1]*m2[2][1]-2.0*m2[2][0]*m[2][1]*m[1][2]-
//                    m[1][0]*m2[2][1]*m[2][2]+m[1][2]*m[2][0]*m2[2][1]-
//                    m[1][1]*m[2][0]*m2[2][2]+m[2][1]*m[1][0]*m2[2][2])/den;
//
//  ans[2][1] = s[2]*-(m[0][1]*m[2][2]*m2[2][0]-m[0][2]*m2[2][0]*m[2][1]-
//                     2.0*m2[2][1]*m[0][0]*m[2][2]+m[1][2]*m2[2][2]*m[0][0]+
//                     2.0*m2[2][1]*m[0][2]*m[2][0]+m[1][0]*m2[2][0]*m[2][2]-
//                     m[1][0]*m[0][2]*m2[2][2]-m[1][2]*m[2][0]*m2[2][0]+
//                     m[0][0]*m2[2][2]*m[2][1]-m2[2][2]*m[0][1]*m[2][0])/den;
//
//  ans[2][2] = s[2]*(m[0][1]*m[1][2]*m2[2][0]-m[0][2]*m2[2][0]*m[1][1]-
//                    m[0][0]*m[1][2]*m2[2][1]+m[1][0]*m[0][2]*m2[2][1]-
//                    m[1][1]*m[2][0]*m2[2][0]-m[2][1]*m2[2][1]*m[0][0]+
//                    2.0*m[1][1]*m2[2][2]*m[0][0]+m[2][1]*m[1][0]*m2[2][0]+
//                    m[2][0]*m[0][1]*m2[2][1]-2.0*m2[2][2]*m[1][0]*m[0][1])/den;
//}

//*********************************************************************
// This is the old version of the zewdie potential where the forces and torques
// were determined completely using mathematica, I'm leaving it in here in case
// we ever want to compare the output of the implementations
// - Tests were performed with the new implementation (in zewdie_analytic) and 
//   the forces and torques matched within a tol=1e-9
// Nonetheless, this code is here for reference
//*********************************************************************
double PairZewdie::zewdie_analytic_old(const int i,const int j, double *quatI, 
				       double *quatP1, double *r12, const double rsq, 
                       double *fforce, double *ttor, double *rtor)
{
  double tempv[3], tempv2[3];
  double temp[3][3];
  double temp1,temp2,temp3;
  double quatItmp[4], quatP1tmp[4];

  int *type = atom->type;
  int newton_pair = force->newton_pair;
  int nlocal = atom->nlocal;


  //josh 
  double potl;
  double uinit1[3],u1[3];
  double uinit2[3],u2[3];
  //assume that the orientation of the particles 
  uinit1[0] = 1; uinit1[1] = 0; uinit1[2] = 0; 
  uinit2[0] = 1; uinit2[1] = 0; uinit2[2] = 0;

  quatItmp[0] = quatI[0]; quatItmp[1] = quatI[1]; quatItmp[2] = quatI[2]; quatItmp[3] = quatI[3];
  quatP1tmp[0] = quatP1[0]; quatP1tmp[1] = quatP1[1]; quatP1tmp[2] = quatP1[2]; quatP1tmp[3] = quatP1[3];
  
  LAMMPS_NS::quat_vec_rot(u1,uinit1,quatItmp);
  LAMMPS_NS::vec_norm(u1);
  LAMMPS_NS::quat_vec_rot(u2,uinit2,quatP1tmp);
  LAMMPS_NS::vec_norm(u2);

  //copy my values into values for andres' code
  double rx,ry,rz, u1x, u1y, u1z, u2x, u2y, u2z;
  rx = r12[0]; ry = r12[1]; rz = r12[2];
  //u1x = u1[0]; u1y = u1[1]; u1z = u1[2]; //define u in for loop, u1x = u1 then u1x = u2
  //u2x = u2[0]; u2y = u2[1]; u2z = u2[2];

  //Define parameters for Zewdie model
  double ps0, ps000, pscc2, ps220, ps222, ps224;
  double pe0, pe000, pecc2, pe220, pe222, pe224;
  //now set them to values from Rippe
  ps0 = 1.0; ps000 = 1.6957; pscc2 = -0.7641; ps220 = -0.1480; ps222 = -0.2582; ps224 = 0.5112;
  pe0 = 1.0; pe000 = 1.9152; pecc2 =  2.7322; pe220 =  1.2633; pe222 =  2.3440; pe224 = 1.0101;
  //these are the values andres used
  //ps0 = 1.0; ps000 = 2.34  ; pscc2 = -1.52  ; ps220 = -0.64  ; ps222 = -0.69  ; ps224 = 1.97  ;
  //pe0 = 1.0; pe000 = 2.24  ; pecc2 =  3.58  ; pe220 =  3.16  ; pe222 =  4.29  ; pe224 = 1.30  ;
 
//  TEST LOOP  
//  //lets calculate the potential as a function of r for diffrent orientations
//  // this code is to check that I've implemented everything correctly
//  int iorient;
//  double rad;
//  ps0 = 1.0;
//  for (iorient=0;iorient < 3; iorient++){
//      if (iorient == 0){
//        u1[0] = 1; u1[1] = 0; u1[2] = 0;
//        u2[0] = 1; u2[1] = 0; u2[2] = 0;
//      }
//      else if (iorient == 1){
//        u1[0] = 0; u1[1] = 0; u1[2] = 1;
//        u2[0] = 0; u2[1] = 0; u2[2] = 1;
//      }
//      else if (iorient == 2){
//        u1[0] = 1; u1[1] = 0; u1[2] = 0;
//        u2[0] = 0; u2[1] = 0; u2[2] = 1;
//      }
//
//      FILE *fptr;
//      char fnme[512];
//      sprintf(fnme,"myu.%d.dat",iorient);
//      fptr = fopen(fnme,"w");
//
//      for(rad=1; rad<5;rad += 0.05){
//        rx = rad; ry = 0.0; rz = 0.0;
//  TEST LOOP  

 
  int ii;
  double norm;
  double drx, dry, drz, dux, duy, duz;
  double Arx, Brx, Crx, Drx, Erx, Zrx, Yrx, Xrx, Wrx, Vrx, Urx, Srx, Rrx, Prx;
  double Ary, Bry, Cry, Dry, Ery, Zry, Yry, Xry, Wry, Vry, Ury, Sry, Rry, Pry;
  double Arz, Brz, Crz, Drz, Erz, Zrz, Yrz, Xrz, Wrz, Vrz, Urz, Srz, Rrz, Prz;
  double Aux, Bux, Cux, Dux, Eux, Zux, Yux, Xux, Wux, Vux, Uux, Sux, Rux, Pux;
  double Auy, Buy, Cuy, Duy, Euy, Zuy, Yuy, Xuy, Wuy, Vuy, Uuy, Suy, Ruy, Puy;
  double Auz, Buz, Cuz, Duz, Euz, Zuz, Yuz, Xuz, Wuz, Vuz, Uuz, Suz, Ruz, Puz;
  double crossx, crossy, crossz;

//for loop to calculate rtor and ttor
  for (ii = 0; ii <2; ii++){
    if (ii==0){
      u1x = u1[0]; u1y = u1[1]; u1z = u1[2];
      u2x = u2[0]; u2y = u2[1]; u2z = u2[2];
    }
    else if (ii==1){
      u1x = u2[0]; u1y = u2[1]; u1z = u2[2];
      u2x = u1[0]; u2y = u1[1]; u2z = u1[2];
    }
    
    //the rx, ry and rz terms are symmetric in u1 and u2, only need to calculate them once...I think
    if (0==ii){
    //code block from andres
    norm = sqrt(pow(rx,2) + pow(ry,2) + pow(rz,2));
    Arx = (rx*u1x + ry*u1y + rz*u1z)/norm;
    Brx = (rx*u2x + ry*u2y + rz*u2z)/norm;
    Crx = u1x*u2x + u1y*u2y + u1z*u2z;
    Drx = (pow(norm,2)*u1x - rx*(rx*u1x + ry*u1y + rz*u1z))/pow(norm,3);
    Erx = (pow(norm,2)*u2x - rx*(rx*u2x + ry*u2y + rz*u2z))/pow(norm,3);
    Zrx = -6*Arx*Drx + 9*Brx*Crx*Drx - 6*Brx*Erx + 9*Arx*Crx*Erx;
    Yrx = 2 - 3*pow(Arx,2) - 3*pow(Brx,2) + 9*Arx*Brx*Crx - 3*pow(Crx,2);
    Xrx = 1 - 5*pow(Brx,2) + 5*pow(Arx,2)*(-1 + 7*pow(Brx,2)) - 20*Arx*Brx*Crx + 2*pow(Crx,2);
    Wrx = 10*(7*pow(Arx,2)*Brx*Erx - Brx*(2*Crx*Drx + Erx) + Arx*((-1 + 7*pow(Brx,2))*Drx - 2*Crx*Erx));
    Vrx = (-2 + 3*pow(Arx,2) + 3*pow(Brx,2))/(2.*sqrt(5));
    Urx = (3*(Arx*Drx + Brx*Erx))/sqrt(5);
    Srx = pscc2*Urx + (ps224*Wrx + 4*ps222*Zrx)/(4.*sqrt(70));
    Rrx = ps000 + ((-1 + 3*pow(Crx,2))*ps220)/(2.*sqrt(5)) + pscc2*Vrx + (ps224*Xrx)/(4.*sqrt(70)) + (ps222*Yrx)/sqrt(70);
    Prx = norm + ps0 - ps0*Rrx;
    drx = (pe0*pow(ps0,6)*((6*(pow(Prx,6) - 2*pow(ps0,6))*(rx - norm*ps0*Srx)*(280*pe000 + 28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx)))/norm + Prx*(-pow(Prx,6) + pow(ps0,6))*(280*pecc2*Urx + sqrt(70)*(pe224*Wrx + 4*pe222*Zrx))))/(70.*pow(Prx,13));
    Ary = -((rx*ry*u1x - pow(norm,2)*u1y + pow(ry,2)*u1y + ry*rz*u1z)/pow(norm,3));
    Bry = -((rx*ry*u2x - pow(norm,2)*u2y + pow(ry,2)*u2y + ry*rz*u2z)/pow(norm,3));
    Cry = norm + ps0 - ps0*Rrx;
    Zry = -6*Arx*Ary - 6*Brx*Bry + 9*Ary*Brx*Crx + 9*Arx*Bry*Crx;
    Yry = 10*(7*pow(Arx,2)*Brx*Bry - Brx*(Bry + 2*Ary*Crx) + Arx*(Ary*(-1 + 7*pow(Brx,2)) - 2*Bry*Crx));
    Xry = (3*(Arx*Ary + Brx*Bry))/sqrt(5);
    Wry = pscc2*Xry + (ps224*Yry + 4*ps222*Zry)/(4.*sqrt(70));
    dry = (pe0*pow(ps0,6)*((6*(pow(Prx,6) - 2*pow(ps0,6))*(ry - norm*ps0*Wry)*(280*pe000 + 28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx)))/norm + Prx*(-pow(Prx,6) + pow(ps0,6))*(280*pecc2*Xry + sqrt(70)*(pe224*Yry + 4*pe222*Zry))))/(70.*pow(Prx,13));
    Arz = -((rx*rz*u1x + ry*rz*u1y - pow(norm,2)*u1z + pow(rz,2)*u1z)/pow(norm,3));
    Brz = -((rx*rz*u2x + ry*rz*u2y - pow(norm,2)*u2z + pow(rz,2)*u2z)/pow(norm,3));
    Zrz = -6*Arx*Arz - 6*Brx*Brz + 9*Arz*Brx*Crx + 9*Arx*Brz*Crx;
    Yrz = 10*(7*pow(Arx,2)*Brx*Brz - Brx*(Brz + 2*Arz*Crx) + Arx*(Arz*(-1 + 7*pow(Brx,2)) - 2*Brz*Crx));
    Xrz = (3*(Arx*Arz + Brx*Brz))/sqrt(5);
    Wrz = pscc2*Xrz + (ps224*Yrz + 4*ps222*Zrz)/(4.*sqrt(70));
    drz = (pe0*pow(ps0,6)*((6*(pow(Prx,6) - 2*pow(ps0,6))*(rz - norm*ps0*Wrz)*(280*pe000 + 28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx)))/norm + Prx*(-pow(Prx,6) + pow(ps0,6))*(280*pecc2*Xrz + sqrt(70)*(pe224*Yrz + 4*pe222*Zrz))))/(70.*pow(Prx,13));
    }

    Aux = (-6*Arx*rx + 9*Brx*Crx*rx + 9*Arx*Brx*norm*u2x - 6*Crx*norm*u2x)/norm;
    Bux = (4*Crx*(-5*Brx*rx + norm*u2x) + 10*Arx*((-1 + 7*pow(Brx,2))*rx - 2*Brx*norm*u2x))/norm;
    Zux = (4*sqrt(14)*Aux*pe222 + sqrt(14)*Bux*pe224 + (168*(Arx*pecc2*rx + Crx*norm*pe220*u2x))/norm)/(56.*sqrt(5));
    Yux = (4*sqrt(14)*Aux*ps222 + sqrt(14)*Bux*ps224 + (168*(Arx*pscc2*rx + Crx*norm*ps220*u2x))/norm)/(56.*sqrt(5));
    dux = (4*pe0*pow(ps0,6)*(-6*ps0*(pow(Prx,6) - 2*pow(ps0,6))*(pe000 + (28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx))/280.)*Yux + Prx*(-pow(Prx,6) + pow(ps0,6))*Zux))/pow(Prx,13);
    Auy = (-6*Arx*ry + 9*Brx*Crx*ry + 9*Arx*Brx*norm*u2y - 6*Crx*norm*u2y)/norm;
    Buy = (4*Crx*(-5*Brx*ry + norm*u2y) + 10*Arx*((-1 + 7*pow(Brx,2))*ry - 2*Brx*norm*u2y))/norm;
    Zuy = (4*sqrt(14)*Auy*pe222 + sqrt(14)*Buy*pe224 + (168*(Arx*pecc2*ry + Crx*norm*pe220*u2y))/norm)/(56.*sqrt(5));
    Yuy = (4*sqrt(14)*Auy*ps222 + sqrt(14)*Buy*ps224 + (168*(Arx*pscc2*ry + Crx*norm*ps220*u2y))/norm)/(56.*sqrt(5));
    duy = (4*pe0*pow(ps0,6)*(-6*ps0*(pow(Prx,6) - 2*pow(ps0,6))*(pe000 + (28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx))/280.)*Yuy + Prx*(-pow(Prx,6) + pow(ps0,6))*Zuy))/pow(Prx,13);
    Auz = (-6*Arx*rz + 9*Brx*Crx*rz + 9*Arx*Brx*norm*u2z - 6*Crx*norm*u2z)/norm;
    Buz = (4*Crx*(-5*Brx*rz + norm*u2z) + 10*Arx*((-1 + 7*pow(Brx,2))*rz - 2*Brx*norm*u2z))/norm;
    Zuz = (4*sqrt(14)*Auz*pe222 + sqrt(14)*Buz*pe224 + (168*(Arx*pecc2*rz + Crx*norm*pe220*u2z))/norm)/(56.*sqrt(5));
    Yuz = (4*sqrt(14)*Auz*ps222 + sqrt(14)*Buz*ps224 + (168*(Arx*pscc2*rz + Crx*norm*ps220*u2z))/norm)/(56.*sqrt(5));
    duz = (4*pe0*pow(ps0,6)*(-6*ps0*(pow(Prx,6) - 2*pow(ps0,6))*(pe000 + (28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx))/280.)*Yuz + Prx*(-pow(Prx,6) + pow(ps0,6))*Zuz))/pow(Prx,13);

    //u1 x du
    crossx = -(u1y*duz - u1z*duy);
    crossy = -(u1z*dux - u1x*duz);
    crossz = -(u1x*duy - u1y*dux);

    if (ii==0){

      ttor[0] = crossx;
      ttor[1] = crossy;
      ttor[2] = crossz;
    }
    else if (ii==1){
      rtor[0] = crossx;
      rtor[1] = crossy;
      rtor[2] = crossz;
    }
  }
  //Now calculate potential
  //potl = (pe0*pow(ps0,6)*(-pow(Prx,6) + pow(ps0,6))*(280*pe000 + 28*sqrt(5)*(-1 + 3*pow(Crx,2))*pe220 + 280*pecc2*Vrx + sqrt(70)*(pe224*Xrx + 4*pe222*Yrx)))/(70.*pow(Prx,12));
  potl =4*pe0*(-(pow(ps0,6)/pow(Prx,6)) + pow(ps0,12)/pow(Prx,12))*(pe000 + ((-1 + 3*pow(Crx,2))*pe220)/(2.*sqrt(5)) + pecc2*Vrx + (pe224*Xrx)/(4.*sqrt(70)) + (pe222*Yrx)/sqrt(70));
  //end code block
  

  //now copy calculated values into returned vectors
  fforce[0] = drx;
  fforce[1] = dry;
  fforce[2] = drz; 

//  TEST LOOP  
//        fprintf(fptr,"%f %f\n",rad,potl);
//      }
//      fclose(fptr);
//  }
//  exit(1);
//  TEST LOOP  

  return potl;
 //end josh
}
