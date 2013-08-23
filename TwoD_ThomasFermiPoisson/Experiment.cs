using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    enum Density_Method
    {
        by_k,
        by_E,
        by_ThomasFermi
    }

    class Experiment
    {
        bool converged = false;
        double alpha, alpha_prime, tol;
        int potential_mixing_rate = 10;

        bool using_flexPDE = false;
        string flexdPDE_input;

        SpinResolved_DoubleMatrix density;
        Band_Data potential;
        DoubleVector band_structure;

        double temperature = 100.0;

        double dy, dz;
        int ny, nz; 

        Density_Method calculation_method = Density_Method.by_ThomasFermi;

        //double fermi_Energy, init_Energy, no_kB_T;
        //double max_Energy;

        //double dk;

        //double dE;
        //int no_Energy_Steps;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dy")) this.dy = (double)input_dict["dy"]; else throw new KeyNotFoundException("No dy in input dictionary!");
            if (input_dict.ContainsKey("ny")) this.ny = (int)(double)input_dict["ny"]; else throw new KeyNotFoundException("No ny in input dictionary!");
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

            // solver inputs
            if (input_dict.ContainsKey("tolerance")) this.tol = (double)input_dict["tolerance"]; else throw new KeyNotFoundException("No solution tolerance in input dictionary!");
            if (input_dict.ContainsKey("alpha")) { this.alpha = (double)input_dict["alpha"]; alpha_prime = alpha; } else throw new KeyNotFoundException("No potential mixing parameter, alpha, in input dictionary!");
            if (input_dict.ContainsKey("FlexPDE_file")) { this.using_flexPDE = true; this.flexdPDE_input = (string)input_dict["FlexPDE_file"]; }

            // physical inputs
            if (input_dict.ContainsKey("T")) this.temperature = (double)input_dict["T"]; else throw new KeyNotFoundException("No temperature in input dictionary!");

            // get the band structure
            if (input_dict.ContainsKey("BandStructure_File")) band_structure = Input_Band_Structure.GetBandStructure((string)input_dict["BandStructure_File"], nz, dz);
            else throw new KeyNotFoundException("No band structure file found in input dictionary!");

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.density = (SpinResolved_DoubleMatrix)input_dict["SpinResolved_Density"]; else this.density = new SpinResolved_DoubleMatrix(ny, nz);
            if (input_dict.ContainsKey("Potential")) this.potential = new Band_Data((DoubleMatrix)input_dict["Potential"]); else potential = Input_Band_Structure.Expand_BandStructure(band_structure / 2.0, ny);

            /*nx = (int)(double)input_dict["nx"];
            dx = (double)input_dict["dx"];

            fermi_Energy = (double)input_dict["E_f"];
            no_kB_T = (double)input_dict["No_kB_T_Above_E_f"];
            temperature = (double)input_dict["T"];


            // momentum density solver parameters
            dk = (double)input_dict["dk"];

            // energy density solver parameters
            init_Energy = (double)input_dict["Lowest_Energy"];
            dE = (double)input_dict["dE"];
            no_Energy_Steps = (int)((max_Energy - init_Energy) / dE);*/
        }

        public void Run()
        {
            DoubleVector donors = new DoubleVector(nz);
            // and put in some delta-dopants
            for (int k = 8; k < 9; k++)
                donors[k] = 0.001;


            // create classes and initialise
            TwoD_PoissonSolver pois_solv = new TwoD_PoissonSolver(dy, dz, ny, nz, using_flexPDE, flexdPDE_input, tol);
            TwoD_ThomasFermiSolver dens_solv = new TwoD_ThomasFermiSolver(band_structure, new DoubleVector(nz), donors, new DoubleVector(nz, -1000.0), new DoubleVector(nz, 1000.0), 0.0, temperature, 10.0, dy, dz, ny, nz);
            
            // run self consistent loop
            int count = 0;
            while (!converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString());

                // calculate the total density for this potential
                density = dens_solv.Get_OneD_Density(potential.mat);

                // solve the potential for the given density and mix in with the old potential
                Band_Data new_potential = (1.0 - alpha) * potential + new Band_Data(alpha * pois_solv.Get_Potential(density.Spin_Summed_Matrix));

                // check for convergence
                converged = pois_solv.Check_Convergence(potential, new_potential, tol);

                // change the potential mixing parameter
                if ((count + 1) % potential_mixing_rate == 0)
                    alpha = pois_solv.Renew_Mixing_Parameter(potential, new_potential, alpha_prime, alpha);

                // transfer new potential array to potential array
                potential = new_potential;
                count++;

                /*/ Calculate 2D electrostatic potential
                TwoD_PoissonSolver pot_solv = new TwoD_PoissonSolver();
                DoubleVector well_potential = pot_solv.Get_Well_Potential(nx);

                // Calculate density using given method
                TwoD_DensitySolver dens_solv = new TwoD_DensitySolver(dx, fermi_Energy, nx);
                dens_solv.Initialise_Hamiltonian(well_potential);

                if (calculation_method == Density_Method.by_k)
                    dens_solv.Solve_Density_Using_Momentum(fermi_Energy, dk);
                else if (calculation_method == Density_Method.by_E)
                    dens_solv.Solve_Density_Using_Energy(dE, init_Energy, no_Energy_Steps);
                else
                    throw new Exception("Error - Method for solving the density not implemented!");

                // calculate DFT potential

                // check if converged
                converged = true;*/
            }

        }


    }
}
