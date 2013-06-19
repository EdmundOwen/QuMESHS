using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace OneD_ThomasFermiPoisson
{
    class Experiment
    {
        bool converged = false;
        double alpha, alpha_prime, tol;
        int potential_mixing_rate;

        double temperature = 100.0;

        double dz;
        int nz;

        public void Initialise(double dz, double alpha, double tol, int nz)
        {
            this.alpha = alpha; this.alpha_prime = alpha;
            this.potential_mixing_rate = 100;

            this.tol = tol;

            this.nz = nz; this.dz = dz;
        }

        public void Run()
        {
            // temporary band structure... (spin-degenerate)
            DoubleVector band_structure = new DoubleVector(nz, 2100.0);
            for (int k = 15; k < nz; k++)
            {
            //int k = 10;
                band_structure[k] = 1400.0;
            }

            DoubleVector donors = new DoubleVector(nz);
            // and put in some delta-dopants
            //for (int k = 15; k < 10; k++)
            //k = 2;
                donors[10] = 0.001;


            // create classes and initialise
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, 0.0, 0.0);
            OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(band_structure, new DoubleVector(nz), donors, new DoubleVector(nz, -1000.0), new DoubleVector(nz, 1000.0), 0.0, temperature, 1.0, nz);

            DoubleVector density = new DoubleVector(2 * nz, 0.0);
            DoubleVector potential = new DoubleVector(nz);
            for (int i = 0; i < nz; i++)
                potential[i] = band_structure[i] / 2.0;

            int count = 0;
            while (!converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString());

                // calculate the total density for this potential
                density = dens_solv.Get_OneD_Density(potential);

                // solve the potential for the given density and mix in with the old potential
                DoubleVector new_potential = potential + alpha * pois_solv.Get_Potential(density);

                // check for convergence
                converged = Check_Convergence(potential, new_potential);

                // change the potential mixing parameter
                if ((count + 1) % potential_mixing_rate == 0)
                    alpha = Renew_Mixing_Factor(potential, new_potential);

                // transfer new potential array to potential array
                potential = new_potential;
                count++;
            }
        }

        /// <summary>
        /// Checks whether the potential has converged by comparing old and new potentials
        /// and determining whether every term is the same within a given tolerance
        /// </summary>
        bool Check_Convergence(DoubleVector potential, DoubleVector new_potential)
        {
            DoubleVector pot_diff = new_potential - potential;

            int[] converged_test = new int[pot_diff.Length];
            for (int i = 0; i < nz; i++)
            {
                converged_test[i] = 0;
                if (Math.Abs(pot_diff[i]) < tol)
                    converged_test[i] = 1;
            }

            if (converged_test.Sum() == pot_diff.Length)
                return true;
            else
                return false;
        }

        double Renew_Mixing_Factor(DoubleVector potential, DoubleVector new_potential)
        {
            // fill a vector with the absolute differences between potential and new_potential (except boundaries)
            DoubleVector pot_diff = new DoubleVector(nz - 2);
            for (int i = 1; i < nz - 1; i++)
                pot_diff[i - 1] = Math.Abs(potential[i] - new_potential[i]);

            // the new mixing factor is the maximum absolute change, multiplied by the original mixing factor
            double new_mixing_parameter = alpha_prime * pot_diff.Max();
            if (new_mixing_parameter < alpha_prime)
                return new_mixing_parameter;
            else
                return alpha_prime;
        }

    }
}
