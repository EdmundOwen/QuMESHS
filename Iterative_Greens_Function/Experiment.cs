using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenCLNet;

namespace Iterative_Greens_Function
{
    class Experiment
    {
        bool converged = false;
        int no_Energy_Count, no_Slices;
        double dE;

        Mem Potential;

        public void Run()
        {
            // Initialise the Experiment

            // Iterate until converged
            while (!converged)
            {
                // Initialise Green's function iterator
                Greens_Function_Calc iter = new Greens_Function_Calc();
                iter.Initialise();

                // Calculate Green's functions
                Mem[] G = iter.Iterate();

                // Calculate Density

                // Solve Poisson's Equation

                // Calculate new potential

                // Check convergence
            }
        }

        /// <summary>
        /// Calculates the density using the diagonal elements of the Green's functions
        /// </summary>
        Mem Calculate_Density(Mem[] G)
        {
        }

        /// <summary>
        /// Calculates DFT potential
        /// </summary>
        Mem Calculate_Potential(Mem Poisson_potential, Mem Density)
        {
        }
    }
}
