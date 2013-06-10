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
        double dE, init_Energy, dx;

        public void Initialise_Experiment(Dictionary<string, object> input_dict)
        {
            no_Slices = (int)(double)input_dict["no_Slices"];
            slice_width = (int)(double)input_dict["slice_width"];
            dx = (double)input_dict["dx"];

            init_Energy = (double)input_dict["initial_Energy"];
            no_Energy_Count = (int)(double)input_dict["no_Energy_Count"];
            dE = (double)input_dict["dE"];
        }

        public void Run()
        {
            // Initialise the Experiment

            Greens_Function_Calc iter = new Greens_Function_Calc(dx, no_Slices, slice_width);
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
                    double Energy = init_Energy + i * dE;
                    Console.WriteLine("Energy =\t" + Energy.ToString());

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
            DoubleMatrix result = new DoubleMatrix(slice_width, no_Slices, 0.0);

            //for (int i = 0; i < no_Slices; i++)
            //{
            //    result[0, i] = 100000;
            //    result[slice_width - 1, i] = 100000;
            //}

            return result;
        }

        /// <summary>
        /// Calculates the density using the diagonal elements of the Green's functions
        /// Warning! INPUT IS OVERWRITTEN
        /// </summary>
        DoubleMatrix Calculate_Density(DoubleComplexMatrix[] G)
        {
            DoubleMatrix result = new DoubleMatrix(slice_width, no_Slices, 0.0);

            for (int i = 0; i < slice_width; i++)
                for (int j = 0; j < no_Slices; j++)
                    result[i, j] = G[j][i, i].Imag / Math.PI;

            return result;
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
