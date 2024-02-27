# Hi-C Metainference using plumed and lammps

## Description

This folder briefly describes how to run a polymer simulation with Hi-C restraints using the metainference method implemented into PLUMED.  
To do this, we have to add to PLUMED an extra cpp files that implement the Hi-C collective variable, and patch the latest version of LAMMPS with PLUMED.

## Installation

- Download PLUMED-2.7.0 from the website  
  <https://github.com/plumed/plumed2/releases/tag/v2.7.0>
- Add the files found in the folder "patches" to the plumed/src/isdb folder ("MetainferenceBase.h", "MetainferenceBase.cpp" and "Metainference.cpp" replace existing files), and remove the files "CS2Backbone.cpp", "EMMI.cpp", "FretEfficiency.cpp", "Jcoupling.cpp", "NOE.cpp", "PRE.cpp", "RCD.cpp", and "SAXS.cpp" from the same directory (the usage of these collective variables have been compromised after modifying metainference to efficiently work with Hi-C data).
- The patches "HIC.cpp" and "HIC2.cpp" implement the contacts to compute the Hi-C map from the simulation using two different kernels. The first is activated with the plumed keyword "HIC" and uses the same kernel introduced in (Carstens, Nilges, and Habeck. "Bayesian inference of chromatin structure ensembles from population-averaged contact data." PNAS 117.14 (2020): 7824). For example:
  ```
  HIC ...
  GAMMA=30
  DC=0.14
  ATOMS1=1,2 REFERENCE1=250.0
  ATOMS2=1,3 REFERENCE2=90.0
  ATOMS3=2,3 REFERENCE3=150.0
  LABEL=hic
  ... HIC
  ```
  defines a Hi-C map with contacts between atoms 1, 2 and 3, with 250, 90, and 150 reads. To compute the contacts, we use the kernel function defined in the paper with a parameter gamma=30 nm$^{-1}$ and dc=0.14 nm. Note that lammps used Angstrom as distant unit by default.  
  On the other hand, the keyword "HIC2" uses a different kernel: $s(d) = \frac{ 1 - ((d-d_0)/r_0)^6 } {1 - ((d-d_0)/r_0)^{12} }$, where $d$ ist the distance between two polymer beads, $d_0$ sets the width of the contact, and $r_0$ sets its slope. This kernel has been historically used by the metadynamics community. For example, the Hi-C map can be defined with:
    ```
  HIC2 ...
  R_0=0.1
  D_0=0.05
  ATOMS1=1,2 REFERENCE1=250.0
  ATOMS2=1,3 REFERENCE2=90.0
  ATOMS3=2,3 REFERENCE3=150.0
  LABEL=hic
  ... HIC2
  ```
  Testing showed that the HIC2 kernel is recommended, since the long-range nature of the PNAS kernel may lead to artifacts.
- Install PLUMED using the installation guidelines. Normally this is simply:  
  ```bash
    ./configure --prefix=/YOUR/PLUMED/INSTALLATION/FOLDER
    source ./sourceme.sh
    make; make install
  ```
   But please check the details as extra flags may be required depending on the machine.
- Remember to let your system know where to find PLUMED. For example, you can do this by adding at the end of your "~/.bashrc" file (if your shell uses bash) the lines:  
  `export PATH="$PATH:/PLUMED/INSTALLATION/FOLDER/bin"`  
  `export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:/PLUMED/INSTALLATION/FOLDER/lib"`  
  `export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/PLUMED/INSTALLATION/FOLDER/lib"`  
  `export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/PLUMED/INSTALLATION/FOLDER/lib/pkgconfig"`  
  `export PLUMED_KERNEL="/PLUMED/INSTALLATION/FOLDER/lib/libplumedKernel.so"`  
- Download LAMMPS-2021 from the website and install with the plumed package active. From within the LAMMPS folder you can do:
  ```bash
    mkdir build-w-plumed
    cd build-w-plumed
    cmake ../cmake -DPKG_USER-PLUMED=yes -DPKG_ASPHERE=yes -DPKG_MOLECULE=yes -DPKG_REPLICA=yes
    make
    make install
  ```
  On Velox, add ` -DCMAKE_C_COMPILER=icc  -DCMAKE_CXX_COMPILER=icpc` during cmake.  
  Both PLUMED and REPLICA packages are definetely required to run polymer metainference simulations. Additional packages may be required depending on the used polymer model.  
  If you encounter errors during the lammps installation, you may need either an older or more recent version of lammps.

## Application to a copolymer with 200 beads

In this application, we first generate an unbiased simulation of a copolymer model (in the folder "run_unbiased", the simulation has been already performed). Then, we define the reference hic map using the script "compute_hic.py" using one of the two following commands depending on the chosen kernel:
```bash
# for HIC kernel
./compute_hic.py -s conf.pdb -f run_unbiased/run_unbiased.xtc -gamma 3.0 -dc 1.4 -o hic_gamma3dc1.4.dat -avg avg_hic_gamma3dc1.4.dat
# for HIC2 kernel
./compute_hic.py -s conf.pdb -f run_unbiased/run_unbiased.xtc -gamma 1 -dc 0.5 -o hic_r1d0.5.dat -avg avg_hic_r1d0.5.dat -rational 1
```
Note that lammps and the scripts use Angstrom as distance units, while plumed uses nm, so you have to be careful about the conversion.

Using the reference Hi-C map, we try to reconstruct the copolymer conformations using metainference and a prior polymer model that does not have any information about structure. In the polymer model defined here, each bead represents 5 nucleosoms. The equilibrium distance between 2 coarse-grained beads is set to 1 Angstrom in lammps, but to get the real size of the polymer we should multiply by a factor of 225. This is very important in real applications. The angle potential was parametrized based on a finer chromatin fiber model that matches the correct sedimentation coefficient of 12-nucleosome fibers.

The scripts "hic_setup.sh" and "hic2_stup.sh" create the reference Hi-C plumed input files for the two possible kernel functions. The parameters (gamma, and dc, or r_0 and d_0) have been set so that the kernel functions have about the same width and slope as a standard Lennard-Jones potential. Parameters that vary very much from these are not recommended.

In the two folders "run_gamma3dc1.4_multi32_noscale" and "run_r1d0.5_multi32_noscale" you find the input files to run the Hi-C metainference simulations with 32 replicas using the two possible kernel functions. Testing showed that 32 replicas as sufficient to reveal subpopulations with a probability as small as 0.05. More replicas will be required to reveal even smaller populations. The file "in.run_multi" is the lammps input file containing the details of the prior polymer model and the name of the plumed input file. The main plumed input file is "plumed_multi.dat" and the only difference between the two folders is the definition of the Hi-C contact map, the metainference settings are the same:
```
METAINFERENCE ...
    ARG=(hic\.hic-.*)
    PARARG=(hic\.exp-.*)
    NOISETYPE=GAUSS
    OPTSIGMAMEAN=SEM AVERAGING=200
    #OPTSIGMAMEAN=NONE # not sampled
    SIGMA_MEAN0=1.0
    SIGMA0=100.0 SIGMA_MIN=0.0001 SIGMA_MAX=1000.0 DSIGMA=0.1
    SCALEDATA SCALE_PRIOR=FLAT SCALE0=300.0 DSCALE=1.0 SCALE_MIN=10.0 SCALE_MAX=1000.0
    WRITE_STRIDE=1000
    LABEL=meta
    TEMP=300.0 # only needed if code does not passes temperature to plumed
... METAINFERENCE

```
The settings mean that we use a single Gaussian prior for the error model (all the Hi-C contacts are assumed to have the same error), the error on the mean Hi-C is computed on the fly during the last 200 steps of the simulation, and both the systematic error and the scale parameter are sampled using a Monte-Carlo scheme. These settings should work in most appllications but you may have to change the ranges of the scale or the sigma depending on the application.
On Velox, a job can be submitted using the script "qsub.sh".