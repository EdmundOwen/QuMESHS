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
        double alpha, tol;

        double dz;
        int nz;

        public void Initialise(double dz, double alpha, double tol, int nz)
        {
            this.alpha = alpha; this.tol = tol;
            this.nz = nz; this.dz = dz;
        }

        public void Run()
        {
            // temporary band structure... (spin-degenerate)
            DoubleVector band_structure = new DoubleVector(2* nz, 2100.0);
            //for (int k = 20; k < 30; k++)
            //{
            int k = 2;
                band_structure[k] = 1400.0;
                band_structure[k + nz] = band_structure[k];
            //}

            DoubleVector donors = new DoubleVector(nz);
            // and put in some delta-dopants
            //for (int k = 5; k < 10; k++)
            k = 5;
                donors[k] = 1.0;


            // create classes and initialise
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, 0.0, 0.0);
            OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(band_structure, new DoubleVector(nz), donors, 0.0, 0.0, 0.0, 100.0, 1.0, 1, 1, nz);

            DoubleVector density = new DoubleVector(2 * nz, 0.0);
            DoubleVector potential = new DoubleVector(2 * nz);
            for (int i = 0; i < 2 * nz; i++)
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
                DoubleVector pot_dens = new_potential - potential;
                int[] converged_test = new int[2 * nz];
                for (int i = 0; i < 2 * nz; i++)
                {
                    converged_test[i] = 0;
                    if (Math.Abs(pot_dens[i]) < tol)
                        converged_test[i] = 1;
                }

                if (converged_test.Sum() == 2 * nz)
                    converged = true;

                // transfer new potential array to potential array
                potential = new_potential;
                count++;
            }
        }
    }
}
