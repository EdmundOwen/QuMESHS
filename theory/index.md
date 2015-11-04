---
title: Theory
---

# Theory

A detailed description of the theory used in QuMESHS can be found
in the [QuMESHS user guide]() but here we give a brief outline of the
methods used.

### Poisson`s Equation with Finite-Element Analysis

Dopents in a semiconductor heterostructure can ionize such that there 
are internal electric fields within the device.  For a given charge 
density distribution, these fields can be calculated by solving the
Poisson equation

*Poisson equation here*

−∇⋅(ϵ(r⃗ )∇ϕ(r⃗ ))=ρ(r⃗ )ϵ0

The finite-element programs used by QuMESHS allow the user to solve 
this equation for an arbitrary device.  The application of voltages 
on the surface gates are included in the boundary conditions of this 
equation.

### Density functional theory

The charge density of the electrons within a semiconductor device is 
given by the solution to the quantum many-body equation

*Schroedinger equation here*
H^Ψ=EΨ

where Ψ is a many-body wave function.  This equation is impossible 
to solve so QuMESHS uses density functional theory to reduce the 
complexity of this problem.  Instead, for a given electrostatic 
potential, we solve the Kohn-Sham equations

*Kohn-Sham equations here*
(−ℏ2∇22m∗+eϕ(r⃗ )+Vxc[n](r⃗ ))Ψj(r⃗ )=EjΨj(r⃗ )

The potential ϕ(r⃗ ) is offset by the band gap of the material 
which varies across a layered semiconductor heterostructure.

### Adaptive Newton step iteration

For Poisson`s equation, we need a density to find a potential whilst 
for DFT we need a potential to find a charge. QuMESHS finds a 
self-consistent solution to this problem using an adaptive Newton step 
method.
