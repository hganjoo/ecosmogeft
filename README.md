# ECOSMOG-EFT v 1.0.0

Himanish Ganjoo - 9 Apr 2026
Himanish.Ganjoo@obspm.fr

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- [![arXiv](https://img.shields.io/badge/arXiv-2007.03042%20-green.svg)](https://arxiv.org/abs/2007.03042) -->

A open-source N-body simulation code for dark-matter only cosmological structure formation for cubic screening in the Effective Field Theory of Dark Energy. The model is taken from [Cusin et al 2018](https://arxiv.org/abs/1712.02782). It is implemented in a modified version of the [ECOSMOG-CVG](https://arxiv.org/abs/2007.03042) (code [here](https://github.com/Christovis/ecosmog-cvg)) which is based on ECOSMOG and [RAMSES](https://github.com/ramses-organisation/ramses). It uses adaptive mesh refinement and adaptive time integration to simulate self-gravitating fluids and is massively parallelizable as it makes use of the MPI communication library.


The code can simulate EFTofDE with Vainshtein screening (eft = True) or linearised EFTs (eft=True and eftlin=True). 


Namelist examples are provided in the `namelist/` folder. 

## Parameters:

### Run params:

- For EFT, set `eft=.true.`
- For linear EFT, set `eft=.true.` and `eftlin=.true.` This version directly solves the modified Poisson equation for gravity without computing the scalar field chi. 

### EFT params:

- `alphaB0`: braiding parameter at a = 1
- `alphaM0`: running of Planck mass at a = 1
- `w0`, `wa`: DE EoS params for CPL parameterisation w(a) = w_0 + (1-a)w_a
- `npre`, `npost`: pre and post smoothing cycles for EFT solver
- `scaling`: 'de' (proportional to Omega_DE / Omega_DE0) or 'a' (proportional to scale factor)
- `nb`, `nm`: exponents for the alpha-scaling.

## Output options:

For EFT switched on, the code will dump phi, grad(phi), psi, grad(psi), chi, grad(chi) at all particle positions as part of the particle file. 
For this option, add the compile time option -DOUTPUT_EXTRADOF_PART to the makefile

**The phi and psi are the _modified_ versions of the two fields, which involve contributions from the extra scalar field chi.**

If eflin=.true., the code will output the same fields, but these have no meaning since the chi is not computed in the linear case. Psi and Chi can be obtained in this case by rescaling phi by mu_chi and mu_psi (see Cusin et al)

**In TESTING**

## Compile time directory-making:

If compiled with -DWITHOUTMKDIR, directories are not created by the code. You will have to create the output directory structure before running the simulation.

### Build/Install
The code has been tested and run on [Infinity](https://infinity-cluster.projet-horizon.fr/) at the Institut Astrophysique de Paris which runs Rocky Linux 8.10, and the [Irene](https://hpc.cea.fr/tgcc-public/en/html/tgcc-public.html) supercomputing cluster of the TGCC.

**Prerequisites:**
* A Fortran compiler
* CMake
* Intel MPI (2018)
* IntelComp (2018)

To build the Fortran BMI bindings from source with cmake, run

```
$ cd ./src/bin
$ make
```

You can quickly test your installation by executing:
```
$ cd bin
$ make
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

Initial condition can be generated with [2LPTic](https://arxiv.org/abs/astro-ph/0606505) from a power spectrum generated with [CAMB](https://github.com/cmbant/CAMB).

The simulation results can be analysed through [Astrild](https://github.com/Christovis/astrild).
