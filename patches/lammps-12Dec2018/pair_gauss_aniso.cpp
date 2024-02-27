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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_gauss_aniso.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "math_vector.h"
#include "atom_vec_ellipsoid.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define SMALL 0.001

/* ---------------------------------------------------------------------- */

PairGaussAniso::PairGaussAniso(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairGaussAniso::~PairGaussAniso()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(r0);
    memory->destroy(theta0);
    memory->destroy(phi0);
    memory->destroy(Ktheta);
    memory->destroy(Kphi);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairGaussAniso::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **tor = atom->torque;

  double ener,fforce[3],r12[3],torqI[3];
  double *quat;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  int myi, myj;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        //always use the higher type as i, this has to due with how angles are defined
        if (itype > jtype){
            myi=i;
            myj=j;
            r12[0] = delx;
            r12[1] = dely;
            r12[2] = delz;
        }
        else if (itype < jtype){
            myi=j; 
            myj=i;
            r12[0] = -delx;
            r12[1] = -dely;
            r12[2] = -delz;
        }
        else{
          error->all(FLERR,"itype cannot equal jtype for pair_lj_aniso!");
        }

        quat = bonus[ellipsoid[myi]].quat;
        ener = aniso_analytic(myi,myj,quat, r12, fforce, torqI);

        fforce[0] *= factor_lj;
        fforce[1] *= factor_lj;
        fforce[2] *= factor_lj;
        torqI[0] *= factor_lj;
        torqI[1] *= factor_lj;
        torqI[2] *= factor_lj;

        f[myi][0] += fforce[0]; 
        f[myi][1] += fforce[1];
        f[myi][2] += fforce[2];
        tor[myi][0] += torqI[0];
        tor[myi][1] += torqI[1];
        tor[myi][2] += torqI[2];

        if (newton_pair || myj < nlocal) {
          f[myj][0] -= fforce[0];
          f[myj][1] -= fforce[1];
          f[myj][2] -= fforce[2];
        }

        if (eflag) {
          evdwl = factor_lj*ener;
        }

        // orig from lj_cut
        //if (evflag) ev_tally(i,j,nlocal,newton_pair,
        //                     evdwl,0.0,fpair,delx,dely,delz);
        // from pair_zewdie
        if (evflag) ev_tally_xyz(myi,myj,nlocal,newton_pair,
                                 evdwl,0.0,fforce[0],fforce[1],fforce[2],
                                 -r12[0],-r12[1],-r12[2]);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairGaussAniso::compute_inner()
{
  error->all(FLERR,"Not implemented! copy from pair_lj_cut.cpp and fix");

}

/* ---------------------------------------------------------------------- */

void PairGaussAniso::compute_middle()
{
  error->all(FLERR,"Not implemented! copy from pair_lj_cut.cpp and fix");
}

/* ---------------------------------------------------------------------- */

void PairGaussAniso::compute_outer(int eflag, int vflag)
{
  error->all(FLERR,"Not implemented! copy from pair_lj_cut.cpp and fix");
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGaussAniso::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(r0,n+1,n+1,"pair:r0");
  memory->create(theta0,n+1,n+1,"pair:theta0");
  memory->create(phi0,n+1,n+1,"pair:phi0");
  memory->create(Ktheta,n+1,n+1,"pair:Ktheta");
  memory->create(Kphi,n+1,n+1,"pair:Kphi");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGaussAniso::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGaussAniso::coeff(int narg, char **arg)
{
  if (narg < 9 || narg > 10)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);
  double r0_one = force->numeric(FLERR,arg[4]);
  double theta0_one = force->numeric(FLERR,arg[5]);
  double phi0_one = force->numeric(FLERR,arg[6]);
  double Ktheta_one = force->numeric(FLERR,arg[7]);
  double Kphi_one = force->numeric(FLERR,arg[8]);

  theta0_one *= MY_PI / 180.0; //convert to radians
  phi0_one *= MY_PI / 180.0; 

  double cut_one = cut_global;
  if (narg == 10) cut_one = force->numeric(FLERR,arg[9]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      r0[i][j] = r0_one;
      theta0[i][j] = theta0_one;
      phi0[i][j] = phi0_one;
      Ktheta[i][j] = Ktheta_one;
      Kphi[i][j] = Kphi_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGaussAniso::init_style()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec) error->all(FLERR,"Pair zewdie requires atom style ellipsoid");

  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this,instance_me);
    else if (respa == 1) {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this,instance_me);

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}


/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */
// removed for compatability with LAMMPS May 2018
//void PairGaussAniso::init_list(int id, NeighList *ptr)
//{
//  if (id == 0) list = ptr;
//  else if (id == 1) listinner = ptr;
//  else if (id == 2) listmiddle = ptr;
//  else if (id == 3) listouter = ptr;
//}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGaussAniso::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  sigma[j][i] = sigma[i][j];
  epsilon[j][i] = epsilon[i][j];
  r0[j][i] = r0[i][j];
  Ktheta[j][i] = Ktheta[i][j];
  Kphi[j][i] = Kphi[i][j];
  theta0[j][i] = theta0[i][j];
  phi0[j][i] = phi0[i][j];

  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  //if (tail_flag) {
  //  int *type = atom->type;
  //  int nlocal = atom->nlocal;

  //  double count[2],all[2];
  //  count[0] = count[1] = 0.0;
  //  for (int k = 0; k < nlocal; k++) {
  //    if (type[k] == i) count[0] += 1.0;
  //    if (type[k] == j) count[1] += 1.0;
  //  }
  //  MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

  //  double sig2 = sigma[i][j]*sigma[i][j];
  //  double sig6 = sig2*sig2*sig2;
  //  double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
  //  double rc6 = rc3*rc3;
  //  double rc9 = rc3*rc6;
  //  etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
  //    sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
  //  ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
  //    sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  //}

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGaussAniso::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&r0[i][j],sizeof(double),1,fp);
        fwrite(&theta0[i][j],sizeof(double),1,fp);
        fwrite(&phi0[i][j],sizeof(double),1,fp);
        fwrite(&Ktheta[i][j],sizeof(double),1,fp);
        fwrite(&Kphi[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGaussAniso::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&r0[i][j],sizeof(double),1,fp);
          fread(&theta0[i][j],sizeof(double),1,fp);
          fread(&phi0[i][j],sizeof(double),1,fp);
          fread(&Ktheta[i][j],sizeof(double),1,fp);
          fread(&Kphi[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&r0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&theta0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&phi0[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&Ktheta[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&Kphi[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGaussAniso::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGaussAniso::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairGaussAniso::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g %g %g\n",i,epsilon[i][i],sigma[i][i],r0[i][i],theta0[i][i],phi0[i][i],Ktheta[i][i],Kphi[i][i]);

}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairGaussAniso::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],r0[i][i],cut[i][j],theta0[i][j],phi0[i][j],Ktheta[i][j],Kphi[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairGaussAniso::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj;
  error->all(FLERR,"Not implemented!");

  return 1;
}

/* ---------------------------------------------------------------------- */

void *PairGaussAniso::extract(const char *str, int &dim)
{
  dim = 7;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"r0") == 0) return (void *) r0;
  if (strcmp(str,"theta0") == 0) return (void *) theta0;
  if (strcmp(str,"phi0") == 0) return (void *) phi0;
  if (strcmp(str,"Ktheta") == 0) return (void *) Ktheta;
  if (strcmp(str,"Kphi") == 0) return (void *) Kphi;
  return NULL;
}

double PairGaussAniso::aniso_analytic(const int i, const int j, double *quatI, 
				        double *r12, double *fforce, double *torqI)
{
    //i is the site whose f and u vectors are used to compute angles
    //j is only used to compute r

    int *type = atom->type;
    int newton_pair = force->newton_pair;
    int nlocal = atom->nlocal;
    int itype, jtype;
    itype = type[i];
    jtype = type[j];

    // Note: to accelerate this code, one could add statements like:
    //       "if (!(newton_pair || j < nlocal)) {}"
    //     to wrap code blocks that are unnecessary of newton_pair = false
    //     since some evaluations are not used.
    //     currently I dont do this, so newton_pair=false will give effectively 
    //     no computational speedup (though inter-processor communication will decrease, 
    //       which could make this a good idea for certain appli


    // get f and u vectors of site I
    double f0[3],u0[3],fI[3],uI[3],qI[4];
    f0[0] = 1; f0[1] = 0; f0[2] = 0; 
    u0[0] = 0; u0[1] = 0; u0[2] = 1; 
    qI[0] = quatI[0]; qI[1] = quatI[1]; qI[2] = quatI[2]; qI[3] = quatI[3];
    LAMMPS_NS::quat_vec_rot(fI,f0,qI);
    LAMMPS_NS::quat_vec_rot(uI,u0,qI);
    LAMMPS_NS::vec_norm(fI);
    LAMMPS_NS::vec_norm(uI);

    // get r
    double r[3],rh[3],rn;
    r[0] = -r12[0]; 
    r[1] = -r12[1]; 
    r[2] = -r12[2];
    rn = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    rh[0] = r[0]/rn;
    rh[1] = r[1]/rn;
    rh[2] = r[2]/rn;


    //----------------------------------------------
    // now calculate angles for modulation
    //----------------------------------------------
    double theta, phi, dtheta, dphi;
    double udotr, fdotr; //dot products
    double mytheta0 = theta0[jtype][itype]; 
    double myphi0 = phi0[jtype][itype];
    
    //theta from u and rh
    udotr = LAMMPS_NS::vec_dot(rh,uI);
    if (udotr > 1) udotr = 1;
    if (udotr <-1) udotr = -1;
    theta = acos(udotr);
    dtheta = theta - mytheta0;
    
    //theta from f and rh
    fdotr = LAMMPS_NS::vec_dot(rh,fI);
    if (fdotr > 1) fdotr = 1;
    if (fdotr <-1) fdotr = -1;
    phi = acos(fdotr);
    dphi = phi - myphi0;

    //----------------------------------------------
    // which range of theta or phi are you in?
    // calculate f_theta and f_phi;
    //----------------------------------------------
    double myKphi, myKtheta;
    int twindow, pwindow;  
    //Ktheta is only upper diagonal, first index must be > than second index, so use jtype first
    myKtheta = Ktheta[jtype][itype]; //2; 
    myKphi = Kphi[jtype][itype];//8;

    if ((dtheta <= MY_PI / myKtheta / 2.0) && (dtheta >= -MY_PI / myKtheta / 2.0))   twindow = 0;
    else if (((dtheta > MY_PI / myKtheta / 2.0) && (dtheta <= MY_PI / myKtheta)) ||
             ((dtheta < -MY_PI / myKtheta / 2.0) && (dtheta >= -MY_PI / myKtheta)))  twindow = 1;
    else if (((dtheta > MY_PI / myKtheta) && (dtheta <= MY_PI)) ||                  
             ((dtheta < -MY_PI / myKtheta) && (dtheta >= -MY_PI)))                 twindow = 2;
    else{ printf("Error! Invalid dtheta (%f)\n",dtheta); exit(1);}

    if ((dphi <= MY_PI / myKphi / 2.0) && (dphi >= -MY_PI / myKphi / 2.0))   pwindow = 0;
    else if (((dphi > MY_PI / myKphi / 2.0) && (dphi <= MY_PI / myKphi)) ||
             ((dphi < -MY_PI / myKphi / 2.0) && (dphi >= -MY_PI / myKphi)))  pwindow = 1;
    else if (((dphi > MY_PI / myKphi) && (dphi <= MY_PI)) ||                  
             ((dphi < -MY_PI / myKphi) && (dphi >= -MY_PI)))               pwindow = 2;
    else{ printf("Error! Invalid dphi (%f)\n",dphi); exit(1);}

    double f_theta, f_phi;
    if      (twindow==0) f_theta = 1;
    else if (twindow==1) f_theta = 1 - cos(myKtheta*dtheta)*cos(myKtheta*dtheta);
    else if (twindow==2) f_theta = 0;

    if      (pwindow==0) f_phi = 1;
    else if (pwindow==1) f_phi = 1 - cos(myKphi*dphi)*cos(myKphi*dphi);
    else if (pwindow==2) f_phi = 0;

    //----------------------------------------------
    // compute LJ energy and forces
    //----------------------------------------------

    double U_gauss, forc_gauss[3];
    double dr,sigma2_inv;

    dr = (rn-r0[itype][jtype]);
    sigma2_inv = 1.0 / sigma[itype][jtype] / sigma[itype][jtype];
    U_gauss = -epsilon[itype][jtype]*exp(-0.5*dr*dr*sigma2_inv);
    
    // force from gaussian
    forc_gauss[0] = - f_theta * f_phi * U_gauss * dr * sigma2_inv * r[0] / rn; 
    forc_gauss[1] = - f_theta * f_phi * U_gauss * dr * sigma2_inv * r[1] / rn;
    forc_gauss[2] = - f_theta * f_phi * U_gauss * dr * sigma2_inv * r[2] / rn;



    //----------------------------------------------
    // now compute forces and torques from f_theta
    //----------------------------------------------
    double forc_theta[3], torq_theta[3];
    if (twindow==0){ 
        torq_theta[0] = 0.0; torq_theta[1] = 0.0; torq_theta[2] = 0.0;
        forc_theta[0] = 0.0; forc_theta[1] = 0.0; forc_theta[2] = 0.0;
    }
    else if (twindow==2){
        forc_theta[0] = 0.0; forc_theta[1] = 0.0; forc_theta[2] = 0.0;
        torq_theta[0] = 0.0; torq_theta[1] = 0.0; torq_theta[2] = 0.0;
    }
    else if (twindow==1){
        double dUdtheta, dftheta_dtheta;
        dftheta_dtheta = 2*myKtheta*cos(myKtheta*dtheta)*sin(myKtheta*dtheta);
        dUdtheta = U_gauss * f_phi * dftheta_dtheta; //NOTE theta AND phi in this expression

        // dtheta/da = dacos(a)/da
        double inverse;
        inverse = sqrt(1-udotr*udotr);
        if (inverse < SMALL) inverse = SMALL;
        
        // dtheta/dr for forces
        double dthetadr[3];
        // dtheta/dr = dtheta/da * da / dr where a = u dot rh
        dthetadr[0] = -1.0/inverse * (uI[0]/rn - udotr*r[0]/rn/rn);
        dthetadr[1] = -1.0/inverse * (uI[1]/rn - udotr*r[1]/rn/rn);
        dthetadr[2] = -1.0/inverse * (uI[2]/rn - udotr*r[2]/rn/rn);

        forc_theta[0] = dUdtheta * dthetadr[0];
        forc_theta[1] = dUdtheta * dthetadr[1];
        forc_theta[2] = dUdtheta * dthetadr[2];

        // dtheta/du for torques
        double dthetadu[3],dUdu[3];
        dthetadu[0] = -rh[0]/inverse;
        dthetadu[1] = -rh[1]/inverse;
        dthetadu[2] = -rh[2]/inverse;
        dUdu[0] = dUdtheta * dthetadu[0];
        dUdu[1] = dUdtheta * dthetadu[1];
        dUdu[2] = dUdtheta * dthetadu[2];
        // u cross dUdu for torque
        torq_theta[0] = -(uI[1]*dUdu[2] - uI[2]*dUdu[1]);
        torq_theta[1] = -(uI[2]*dUdu[0] - uI[0]*dUdu[2]);
        torq_theta[2] = -(uI[0]*dUdu[1] - uI[1]*dUdu[0]);
    }

    //---------------------------------
    //modulate by f_phi
    //---------------------------------
    double forc_phi[3], torq_phi[3];
    if (pwindow==0){
        torq_phi[0] = 0.0; torq_phi[1] = 0.0; torq_phi[2] = 0.0;
        forc_phi[0] = 0.0; forc_phi[1] = 0.0; forc_phi[2] = 0.0;
    }
    else if (pwindow==2){
        forc_phi[0] = 0.0; forc_phi[1] = 0.0; forc_phi[2] = 0.0;
        torq_phi[0] = 0.0; torq_phi[1] = 0.0; torq_phi[2] = 0.0;
    }
    else if (pwindow==1){
        double dUdphi, dfphi_dphi;
        dfphi_dphi = 2*myKphi*cos(myKphi*dphi)*sin(myKphi*dphi);
        dUdphi = U_gauss * f_theta * dfphi_dphi; //NOTE theta AND phi in this expression

        // dphi/da = dacos(a)/da
        double inverse;
        inverse = sqrt(1-fdotr*fdotr);
        if (inverse < SMALL) inverse = SMALL;
        
        // dphi/dr for forces
        double dphidr[3];
        // dphi/dr = dphi/da * da / dr where a = u dot rh
        dphidr[0] = -1.0/inverse * (fI[0]/rn - fdotr*r[0]/rn/rn);
        dphidr[1] = -1.0/inverse * (fI[1]/rn - fdotr*r[1]/rn/rn);
        dphidr[2] = -1.0/inverse * (fI[2]/rn - fdotr*r[2]/rn/rn);

        forc_phi[0] = dUdphi * dphidr[0];
        forc_phi[1] = dUdphi * dphidr[1];
        forc_phi[2] = dUdphi * dphidr[2];
        
        // dphi/du for torques
        double dphidf[3],dUdf[3];
        dphidf[0] = -rh[0]/inverse;
        dphidf[1] = -rh[1]/inverse;
        dphidf[2] = -rh[2]/inverse;
        dUdf[0] = dUdphi * dphidf[0];
        dUdf[1] = dUdphi * dphidf[1];
        dUdf[2] = dUdphi * dphidf[2];
        //cross product for torque
        torq_phi[0] = -(fI[1]*dUdf[2] - fI[2]*dUdf[1]); //make sure fI x dUdf!
        torq_phi[1] = -(fI[2]*dUdf[0] - fI[0]*dUdf[2]);
        torq_phi[2] = -(fI[0]*dUdf[1] - fI[1]*dUdf[0]);
    }

    double U;
    U = f_theta * f_phi * U_gauss;

    fforce[0] = forc_gauss[0] + forc_theta[0] + forc_phi[0];
    fforce[1] = forc_gauss[1] + forc_theta[1] + forc_phi[1];
    fforce[2] = forc_gauss[2] + forc_theta[2] + forc_phi[2];
    
    torqI[0] = torq_theta[0] + torq_phi[0];
    torqI[1] = torq_theta[1] + torq_phi[1];
    torqI[2] = torq_theta[2] + torq_phi[2];
    double ttotal;
    ttotal = torqI[0]*torqI[0] + torqI[1]*torqI[1] + torqI[2]*torqI[2];

    return U;
   
}

