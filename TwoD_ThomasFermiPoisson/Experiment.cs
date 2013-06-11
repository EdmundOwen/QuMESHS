using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace TwoD_ThomasFermiPoisson
{
    enum Density_Method
    {
        by_k,
        by_E
    }

    class Experiment
    {
        // physical constants
        public const double kB = 0.086173324;               // (meV) (K)^-1

        bool converged = false;
        Density_Method calculation_method;

        double fermi_Energy, init_Energy, no_kB_T;
        double max_Energy;

        double dk;

        double dE;
        int no_Energy_Steps;

        double temperature;

        int ny, nz;
        int layer_index;
        double dy, dz;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            ny = (int)(double)input_dict["ny"];
            dy = (double)input_dict["dy"];

            fermi_Energy = (double)input_dict["E_f"];
            no_kB_T = (double)input_dict["No_kB_T_Above_E_f"];
            temperature = (double)input_dict["T"];
            // calculate the maximum energy to calculate to above the fermi surface
            max_Energy = fermi_Energy + (no_kB_T * kB * temperature);


            // momentum density solver parameters
            dk = (double)input_dict["dk"];

            // energy density solver parameters
            init_Energy = (double)input_dict["Lowest_Energy"];
            dE = (double)input_dict["dE"];
            no_Energy_Steps = (int)((max_Energy - init_Energy) / dE);
        }

        public void Run()
        {
            //

            // run self consistent loop
            while (!converged)
            {
                // Calculate 2D electrostatic potential
                TwoD_PoissonSolver pot_solv = new TwoD_PoissonSolver();
                DoubleVector well_potential = pot_solv.Get_Well_Potential(ny);

                // Calculate density using given method
                TwoD_DensitySolver dens_solv = new TwoD_DensitySolver(dy, fermi_Energy, ny);
                dens_solv.Initialise_Hamiltonian(well_potential);

                if (calculation_method == Density_Method.by_k)
                    dens_solv.Solve_Density_Using_Momentum(max_Energy, dk);
                else if (calculation_method == Density_Method.by_E)
                    dens_solv.Solve_Density_Using_Energy(dE, init_Energy, no_Energy_Steps);
                else
                    throw new Exception("Error - Method for solving the density not implemented!");

                // calculate DFT potential

                // check if converged
                converged = true;
            }

        }


    }
}
