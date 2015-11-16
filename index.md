---
title: Home
---

# QuMESH

## What is QuMESHS?

QuMESHS is a numerical code for solving finding the electron densities in 
quantum mesoscopic systems.  Specifically, it is designed to model layered 
semiconductor heterostructure devices with arbitrary surface gate geometries. 
The electron density and electrostatic potential are calculated using 
self-consistent iterative procedure.  The potential is calculated using 
finite-element analysis and the density is (typically) calculated using 
density functional theory.  Simulations in 1, 2 and 3D can be performed with 
varying degrees of approximation for the density.  Dopents and surface 
charges are included.  We intend the solutions of QuMESHS to be accurate 
representations of the actual experimental devices without the need to resort 
to using model potentials.

QuMESHS is written in C#.  The density calculations and the control code is 
performed in the main body of the QuMESHS code.  For 1D calculations, QuMESHS 
can also calculate the electrostatic potential but for 2 and 3D, we use 
external programs to calculate the electrostatic potential using the 
finite-element method.  At the moment, this is done using either 
[FlexPDE](http://www.pdesolutions.com) or the [deal.II](http://dealii.org) 
FEA library but QuMESHS can use any external program to generate these 
potentials.
