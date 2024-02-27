# Hi-C Metainference using PLUMED and LAMMPS

## 1. Description

Here we briefly describe how to run a polymer simulation with Hi-C restraints using the metainference method implemented into PLUMED.
To do this, we have to add to PLUMED an extra cpp files that implement the Hi-C collective variable, and patch LAMMPS with PLUMED.  
Importantly, we optimized metainference to run with thousands of experimental observables (as typical with Hi-C).
The original PLUMED metainference inplementation scaled with N^2, with N being the number of observables, while our version scales with N.
However, note that our changes to the main metainference routines make this incompatible with the use of the original experimental observables implemented (such as chemical shifts and SAXS data). Therefore, when patching plumed it is also necessary to remove these corresponding files.
We implemented Hi-C metainference in two versions of PLUMED: 2.7.0 and 2.5.7 (the 2.7.0 patch should be compatible with all 2.7.\* versions, and similarly for the 2.5 patch, but this is not guaranteed).
This is useful for the compatibility between LAMMPS and PLUMED.
In particular, the 1CPN chromatin model works in a limited number of relatively old LAMMPS versions, such as LAMMPS-12Dec2018, which is not compatible with the most recent PLUMED versions. LAMMPS-12Dec2018 is compatible with PLUMED-2.5.7, which is one of the first versions of plumed to implement metainference.  
To compare Hi-C data with contact probabilities in MD, we implement two kernel functions, corresponding to the collective variables HIC and HIC2.
The first is the same kernel introducted in [Carstens et al. PNAS 117.14 (2020): 7824].
The second corresponds to the rational coordination number already implemented as a collective variable in PLUMED <https://www.plumed.org/doc-v2.8/user-doc/html/_c_o_o_r_d_i_n_a_t_i_o_n.html>.

## 2. Patches and installation

### Plumed

- Download PLUMED-2.7.0 (for using a relatively recent version of LAMMPS, for example from LAMMPS-23Jun2022)  
  <https://github.com/plumed/plumed2/releases/tag/v2.7.0>  
  or PLUMED-2.5.7 (for running Hi-C metainference with the 1CPN model using LAMMPS-12Dec2018)  
  <https://github.com/plumed/plumed2/releases/tag/v2.5.7>  
- Add the files found in the corresponding folder `patches/plumed-2.x` to the plumed/src/isdb folder (`MetainferenceBase.h` replaces existing file),
and remove the files `CS2Backbone.cpp`, `EMMI.cpp`, `FretEfficiency.cpp`, `Jcoupling.cpp`, `NOE.cpp`, `PRE.cpp`, `RCD.cpp`, and `SAXS.cpp` from the same directory.
  The usage of these collective variables have been compromised after modifying metainference to efficiently work with Hi-C data.
- The patches `HIC.cpp` and `HIC2.cpp` implement the contacts to compute the Hi-C map from the simulation using two different kernels. The first is activated with the PLUMED keyword `HIC` and uses the same kernel introduced in by Carstens et al., while the second is activated with the keyword HIC2 and the kernel is the same as the standard plumed coordination number.
- Install PLUMED following the normal installation guidelines on the PLUMED website. Normally this is simply:  
  `$ ./configure --prefix=/YOUR/PLUMED/INSTALLATION/FOLDER`  
  `$ source ./sourceme.sh`  
  `$ make; make install`  
  But please check the details as extra flags may be required depending on the machine.
- Remember to let your system know where to find PLUMED. For example, you can do this by adding at the end of your `~/.bashrc` file (if your shell uses bash) the lines:  
  `export PATH="$PATH:/PLUMED/INSTALLATION/FOLDER/bin"`  
  `export DYLD_LIBRARY_PATH="$DYLD_LIBRARY_PATH:/PLUMED/INSTALLATION/FOLDER/lib"`  
  `export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/PLUMED/INSTALLATION/FOLDER/lib"`  
  `export PKG_CONFIG_PATH="$PKG_CONFIG_PATH:/PLUMED/INSTALLATION/FOLDER/lib/pkgconfig"`  
  `export PLUMED_KERNEL="/PLUMED/INSTALLATION/FOLDER/lib/libplumedKernel.so"`  

### Lammps

- Download the version of LAMMPS of your choice from <https://download.lammps.org/tars/index.html> (e.g., versions 23Jun2022 or 12Dec2018) and install it with the PLUMED package active.
  Typically, from within the LAMMPS folder you would do:  
  `$ mkdir build-w-plumed`  
  `$ cd build-w-plumed`  
  `$ cmake ../cmake -D CMAKE_INSTALL_PREFIX=/your/install/directory -DPKG_USER-PLUMED=yes -DPKG_ASPHERE=yes -DPKG_MOLECULE=yes -DPKG_REPLICA=yes -DPKG_MANYBODY=yes -DPKG_EXTRA-DUMP=yes`  
  `$ make; make install`  
  On some computer clusters, you may need to add certains options during cmake, for example `-DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc` to use intel compilers.  
  Both the PLUMED and the REPLICA packages are required to run polymer metainference simulations. Additional packages may be required depending on the used polymer model.  
  If you encounter errors during the LAMMPS installation, you may need either an older or more recent version of LAMMPS.
- To use the 1CPN model with Hi-C metainference, you need to patch the LAMMPS-12Dec2018 source code to add specific potentials used in 1CPN. The files to add to the src directory of LAMMPS are found in the `patches/lammps-12Dec2018/` directory. The file `atom.cpp` has also been edited to use the 1CPN model. A description of these patches is given in the 1CPN website <https://1cpn-model.readthedocs.io/en/latest/lammps.html>. These patches were originally developed for LAMMPS-16Mar2018, but they also work for LAMMPS-12Dec2018 due to the minor changes between the two versions. The patches cannot be used on more recent versions of LAMMPS.

## 3. Usage

To run Hi-C metainference, we first define the experimental contact maps in the plumed input file. For example:

`HIC ...`  
`GAMMA=30`  
`DC=0.14`  
`ATOMS1=1,2 REFERENCE1=250.0`  
`ATOMS2=1,3 REFERENCE2=90.0`  
`ATOMS3=2,3 REFERENCE3=150.0`  
`LABEL=hic`  
`... HIC`  

defines a Hi-C map with contacts between atoms 1, 2 and 3, with 250, 90, and 150 reads. To compute the contacts, we use the kernel function defined in [Carstens et al. PNAS 117.14 (2020): 7824] with parameters gamma=30/nm and dc=0.14 nm. Note that LAMMPS uses Angstrom as distance unit when the units are set to real.  
On the other hand, the keyword `HIC2` uses a different kernel: $s(d) = \frac{ 1 - ((d-d_0)/r_0)^6 } {1 - ((d-d_0)/r_0)^{12} }$, where $d$ ist the distance between two polymer beads, $d_0$ sets the width of the contact, and $r_0$ sets its slope. This kernel has been historically used by the metadynamics community. For example, the Hi-C map can be defined with:

`HIC2 ...`  
`R_0=0.1`  
`D_0=0.05`  
`ATOMS1=1,2 REFERENCE1=250.0`  
`ATOMS2=1,3 REFERENCE2=90.0`  
`ATOMS3=2,3 REFERENCE3=150.0`  
`LABEL=hic`  
`... HIC2`  

Once the contact map is defined, we can define the metainference calculation. For example:  

`METAINFERENCE ...`  
`    ARG=(hic\.hic-.*)`  
`    PARARG=(hic\.exp-.*)`  
`    NOISETYPE=GAUSS`  
`    OPTSIGMAMEAN=SEM AVERAGING=200`  
`    SIGMA_MEAN0=1.0`  
`    SIGMA0=100.0 SIGMA_MIN=0.0001 SIGMA_MAX=1000.0 DSIGMA=0.1`  
`    SCALEDATA SCALE_PRIOR=FLAT SCALE0=300.0 DSCALE=1.0 SCALE_MIN=10.0 SCALE_MAX=1000.0`  
`    WRITE_STRIDE=1000`  
`    LABEL=meta`  
`    TEMP=300.0 # only needed if code does not pass temperature to plumed`  
`... METAINFERENCE`  

defines a metainference bias using a previously-defined `hic` contact map, with a Gaussian noise model, error on the means estimate on the fly from the last 200 frames, and systematic errors and scale (to convert contact probabilities into reads) sampled by Monte Carlo.  
It uses the same standard metainference syntax. Please check the relevant documentation on the plumed website.

## 4. Application

In the folder `application`, we describe an application to a 30kb locus on chromosome 5 of mouse embryonic stem cells, based on the Micro-C data by Hsieh et al. [Hsieh et al. Mol. cell 78.3 (2020): 539], first using a prior model at 1kb resolutioni (folder `md_cg1kb_30kb`), and then using the 1CPN model at nucleosome resolution [Lequieu et al. J. Chem. Phys. 150.21 (2019): 215102] (folder `md_1cpn_n150nrl200`).
The reference experimental data at 1kb or 200bp resolution are given in the folder `experimental_data`.
In both cases we use the rational kernel (HIC2).
The folders contain the scripts and the input files necessary to run the simulations.

### 1kb model

We first produce an unbiased MD trajectory to generate the initial conformations for the Metainference simulation, and then run the actual Hi-C metainference simulation using 64 replicas.
Note that LAMMPS uses Angstrom as distance units, while PLUMED uses nm, so you have to be careful about the conversion.
The equilibrium distance between 2 coarse-grained beads is set to 1 Angstrom in LAMMPS, but to get the real size of the polymer we should multiply by a factor of 220.
In the `run_hic_multi64` folder, the script `hic_setup.sh` creates the PLUMED input file. The parameters r_0 and d_0 have been set so that the kernel functions have about the same width and slope as a standard Lennard-Jones potential. Parameters that vary very much from these are not recommended.

### 1CPN model

In our prior 1CPN model, the nucleosome repeat length is 200 bp, corresponding to the resolution of the Micro-C data and about the same as the average nucleosome repeat length estimated from chemical mapping experiments.
Our system is therefore made of 150 nucleosomes.
The sequence of the nucleosomal DNA is set to poly-ttaggg (we set this in the 1CPN LAMMPS input files), which have an intermediate positioning between 601 DNA and poly-A tracts.
To compare the Micro-C reads and the chromatin contacts, the HIC2 contact kernel is applied to the distance between the histone core particles of the corresponding genomic regions.
The jupyter notebook `setup_hic.ipynb` sets up the plumed contact map based on the experimental data, using a kernel with a size of 13 nm.

To ease the convergence of the 1CPN metainference simulations, we use the following protocol:  

- We equilibrate the 150-nucleosome chromatin fibers at 700 K for 10^8 MD steps (64 independent runs, folder: `run_ttaggg_T700`).
- We take the output conformations from the previous step and start steered MD simulations where the centers of mass of 5-nucleosome regions (=1 kb) are pulled towards the corresponding positions observed in the 1-kb metainference simulations (folder: `run_fit_cg64_coverage`). Check `setup_fit_cg.ipynb` notebook for setup.
- At last, we use the final configurations in the previous step as starting point for the actual metainference simulations (folder: `run_meta64_coverage13_all`).



