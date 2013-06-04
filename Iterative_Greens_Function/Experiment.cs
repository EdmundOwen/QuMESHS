using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenCLNet;
using CenterSpace.NMath.Core;

namespace Iterative_Greens_Function
{
    class Experiment
    {
        bool converged = false;
        int no_Energy_Count, no_Slices, slice_width;
        double dE;

        DoubleMatrix Potential;

        public void Initialise_Experiment(Dictionary<string, object> input_dict)
        {
            no_Slices = (int)input_dict["no_Slices"];
            slice_width = (int)input_dict["slice_width"];

            no_Energy_Count = (int)input_dict["no_Energy_Count"];
            dE = (double)input_dict["dE"];
        }

        public void Run()
        {
            // Initialise the Experiment

            Greens_Function_Calc iter = new Greens_Function_Calc(no_Slices, slice_width);
            DoubleMatrix Potential = Initialise_Potential();

            // Iterate until converged
            while (!converged)
            {
                // Initialise Green's function iterator
                iter.Initialise();

                DoubleMatrix Density = new DoubleMatrix(slice_width, no_Slices);
                // Iterate for all energies
                for (int i = 0; i < no_Energy_Count; i++)
                {
                    double Energy = i * dE;

                    DoubleComplexMatrix G_0n, G_n0;
                    // Calculate Green's functions
                    DoubleComplexMatrix[] G_ii = iter.Iterate(out G_0n, out G_n0, Potential, Energy);

                    // Calculate Density
                    Density += Calculate_Density(G_ii);
                }

                // Solve Poisson's Equation

                // Calculate new potential

                // Check convergence
            }
        }

        DoubleMatrix Initialise_Potential()
        {
            return new DoubleMatrix(slice_width, no_Slices, 0.0);
        }

        /// <summary>
        /// Calculates the density using the diagonal elements of the Green's functions
        /// </summary>
        DoubleMatrix Calculate_Density(DoubleComplexMatrix[] G)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// Calculates DFT potential
        /// </summary>
        DoubleMatrix Calculate_Potential(DoubleMatrix Poisson_potential, DoubleMatrix Density)
        {
            throw new NotImplementedException();
        }
    }
}
