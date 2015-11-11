# QuMESHS
Quantum Mesoscopic Electronic Semiconductor Heterostructure Solver

## List of Contents

1. Introduction
2. Requirements
  1. Running QuMESHS with NMath
  2. Compiling QuMESHS with NMath
  3. Calculating electrostatic potentials with QuMESHS
    1. Using QuMESHS with FlexPDE
    2. Using QuMESHS with deal.II

## 1. Introduction
The QuMESHS program is designed to calculate self-consistent electrostatic potentials and charge densities for layered semiconductor heterostructures.

All of the information found here should be available on our website at www.qumeshs.org along with documentation and examples.  In this README, comments in this format (TODO ...) are features which will hopefully be added in the future.

## 2.  Requirements
QuMESHS uses external libraries and programs which must be installed and configured for it to compile and execute.  The two main usages of external software is for the linear algebra, which is provided by the NMath library, and the finite element calculations of the electrostatic potential, which is currently configured to run using either FlexPDE or deal.II.

### 2.i. Running QuMESHS with NMath
Currently, QuMESHS calculates wavefunctions, finds zeros and integrates functions using the NMath library provided by Centerspace (http://www.centerspace.net/).  The license key for this library is set during compilation.  Therefore, if you just want to run QuMESHS and not alter the underlying code, you do not need to purchase the NMath library.  You can download the most up-to-date QuMESHS release from the release page on GitHub.  QuMESHS should work when the relevant dll is in the same folder which you are running QuMESHS in.

### 2.ii. Compiling with NMath
This NMath library is proprietary software and compiling executables using this library is not possible without purchasing NMath.  The license key is accessed using the Solver_Bases.License.NMath_License_Key property.  If you want to compile your own copy of QuMESHS, then you need to purchase a license and create a static License class with the NMath_License_Key property.  An example license class can be found in at Solver_Bases.Example_License.  (TODO Remove requirement for using NMath by using an open source library)

### 2.iii. Calculating electrostatic potentials with QuMESHS
In 2D and 3D, QuMESHS calculates the electrostatic potential using finite element analysis.  QuMESHS does not do this itself, but instead calls external programs which save the relevant data to the hard drive which QuMESHS then loads up and uses to calculate self-consistent solutions.  At the moment, electrostatic potentials are calculated using either FlexPDE or deal.II

#### 2.iii.a.  Using QuMESHS with FlexPDE
FlexPDE is proprietary software released by PDE Solutions Inc. (http://www.pdesolutions.com/).  In order to run QuMESHS with FlexPDE, the computer must have FlexPDE installed and a full license must be available.  The accuracy of the student and demo solvers is nowhere near good enough to achieve convergense.  At the moment, QuMESHS creates its own FlexPDE script file using the relevant band structure data and it can only do this for split gate devices.  (TODO Create general script file)

#### 2.iii.b. Using QuMESHS with deal.II
deal.II is open-source software provided by the deal.II team (dealii.org).  This is a finite element library for Linux and must be downloaded and compiled yourself.  At the moment, to use QuMESHS yourself, you must create your own deal.II executable which outputs a potential file with the correct format for QuMESHS to read in.  This is quite a bit of work but the deal.II library is fantastically documented on its website with many examples which make a good framework.  (TODO Add example deal.II files to website) (TODO Create deal.II program for calculating with general device structures).  More details can be found on in the user guide.
