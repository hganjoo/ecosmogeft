# ECOSMOG-EFT

Himanish Ganjoo - 11 Aug 2025

Himanish.Ganjoo@obspm.fr

This code implements the Effective Field Theory of Dark Energy with cubic screening in RAMSES. 

Adapted from the ECOSMOG-CVG code by Christoph Becker (https://github.com/Christovis/ecosmog-cvg). 

EFT theory taken from Cusin et al (https://arxiv.org/abs/1712.02782). 

Namelist examples are provided in the namelist/ folder. 

## Parameters:

### Run params:

- For EFT, set eft=.true.
- For linear EFT, set eft=.true. and eftlin=.true. This version directly solves the modified Poisson equation for gravity without computing the scalar field chi. 

### EFT params:

- alphaB0 = braiding parameter at a = 1
- alphaM0 = running of Planck mass at a = 1
- w0, wa = DE EoS params
- npre, npost = pre and post smoothing cycles for EFT solver

## Output options:

For EFT switched on, the code will dump phi, grad(phi), psi, grad(psi), chi, grad(chi) at all particle positions as part of the particle file. 
For this option, add the compile time option -DOUTPUT_EXTRADOF_PART to the makefile

**The phi and psi are the _modified_ versions of the two fields, which involve contributions from the extra scalar field chi.**

If eflin=.true., the code will output the same fields, but these have no meaning since the chi is not computed in the linear case. Psi and Chi can be obtained in this case by rescaling phi by mu_chi and mu_psi (see Cusin et al)


