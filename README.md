# QuMESHS
Quantum Mesoscopic Electronic Semiconductor Heterostructure Solver

## List of Contents

1. Introduction
2. Requirements
  1. Running QuMESHS with NMath
  2. Compiling QuMESHS with NMath
  3. Using QuMESHS with FlexPDE
  4. Using QuMESHS with deal.II

## 1. Introduction
The QuMESHS program is designed to calculate self-consistent electrostatic 
potentials and charge densities for layered semiconductor heterostructures.

The software was developed for...

## x.  Requirements
QuMESHS uses external libraries and programs which must be installed and
configured for it to compile and execute.  The two main usages of external
software is for the linear algebra, which is provided by the NMath library,
and the finite element calculations of the electrostatic potential, which
is currently configured to run using either FlexPDE or deal.II.

### x.x. Running QuMESHS with NMath
Currently, QuMESHS calculates wavefunctions, finds zeros and integrates
functions using the NMath library provided by Centerspace 
(http://www.centerspace.net/).  The license key for this library is set
during compilation.  Therefore, if you just want to run QuMESHS and not
alter the underlying code, you do not need to purchase the NMath library.
You can download the most up-to-date QuMESHS release from the release
page on GitHub.  QuMESHS should work when the relevant dll is in the same 
folder which you are running QuMESHS in.

### x.x. Compiling with NMath
This NMath library is proprietary software and compiling executables 
using this library is not possible without purchasing NMath.  The license 
key is accessed using the Solver_Bases.License.NMath_License_Key property.  
If you want to compile your own copy of QuMESHS, then you need to purchase 
a license and create a static License class with the NMath_License_Key 
property.  An example license class can be found in at 
Solver_Bases.Example_License.

### x.x. Using QuMESHS with FlexPDE

### x.x. Using QuMESHS with deal.II
