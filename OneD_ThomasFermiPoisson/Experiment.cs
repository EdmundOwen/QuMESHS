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
        int potential_mixing_rate = 100;
        
        bool using_flexPDE = false;
        string flexdPDE_input;

        DoubleVector potential, density, band_structure;

        double temperature = 100.0;

        double dz;
        int nz;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

            // solver inputs
            if (input_dict.ContainsKey("tolerance")) this.tol = (double)input_dict["tolerance"]; else throw new KeyNotFoundException("No solution tolerance in input dictionary!");
            if (input_dict.ContainsKey("alpha")) { this.alpha = (double)input_dict["alpha"]; alpha_prime = alpha; } else throw new KeyNotFoundException("No potential mixing parameter, alpha, in input dictionary!");
            if (input_dict.ContainsKey("FlexPDE_file")) { this.using_flexPDE = true; this.flexdPDE_input = (string)input_dict["FlexPDE_file"]; }

            // physical inputs
            if (input_dict.ContainsKey("T")) this.temperature = (double)input_dict["T"]; else throw new KeyNotFoundException("No temperature in input dictionary!");
            if (input_dict.ContainsKey("BandStructure_File")) band_structure = GetBandStructure((string)input_dict["BandStructure_File"]); else throw new KeyNotFoundException("No band structure file found in input dictionary!");

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.density = (DoubleVector)input_dict["SpinResolved_Density"]; else this.density = new DoubleVector(2 * nz, 0.0);
            if (input_dict.ContainsKey("Potential")) this.potential = (DoubleVector)input_dict["Potential"]; else potential = band_structure / 2.0;
        }

        public void Initialise(DoubleVector Potential, DoubleVector Density, double dz, double alpha, double tol, int nz)
        {
            this.potential = Potential; this.density = Density;
            Initialise(dz, alpha, tol, nz);
        }

        public void Initialise(double dz, double alpha, double tol, int nz)
        {

            this.alpha = alpha; this.alpha_prime = alpha;
            this.potential_mixing_rate = 100;

            this.tol = tol;

            this.nz = nz; this.dz = dz;
        }

        public void Run()
        {

            DoubleVector donors = new DoubleVector(nz);
            // and put in some delta-dopants
            for (int k = 27; k < 43; k++)
                donors[k] = 0.0001;


            // create classes and initialise
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, 0.0, 0.0, using_flexPDE, flexdPDE_input);
            OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(band_structure, new DoubleVector(nz), donors, new DoubleVector(nz, -900.0), new DoubleVector(nz, 900.0), 0.0, temperature, 1.0, dz, nz);
            
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

            Output(density, "density.dat");
            Output(potential, "potential.dat");
        }

        void Output(DoubleVector data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            for (int i = 0; i < data.Length; i++)
                sw.WriteLine(data[i].ToString());

            sw.Close();
        }


        DoubleVector GetBandStructure(string filename)
        {
            // temporary band structure... (spin-degenerate)
            DoubleVector result = new DoubleVector(nz, 2100.0);
            for (int k = 100; k < 130; k++)
            {
                //int k = 10;
                result[k] = 1400.0;
            }
            
            return result;
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
