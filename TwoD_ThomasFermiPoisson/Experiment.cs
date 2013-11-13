using System;
using System.Collections.Generic;
using System.IO;
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

    public class Experiment : Experiment_Base
    {
        SpinResolved_DoubleMatrix charge_density;

        double dy, dz;
        int ny, nz; 

        Density_Method calculation_method = Density_Method.by_ThomasFermi;

        public Experiment()
        {
        }

        public void Initialise_Experiment(Dictionary<string, object> input_dict)
        {
            Console.WriteLine("Initialising Experiment");

            // simulation domain inputs
            if (input_dict.ContainsKey("dy")) this.dy = (double)input_dict["dy"]; else throw new KeyNotFoundException("No dy in input dictionary!");
            if (input_dict.ContainsKey("ny")) this.ny = (int)(double)input_dict["ny"]; else throw new KeyNotFoundException("No ny in input dictionary!");
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");
            
            // physics parameters are done by the base method
            base.Initialise(input_dict, dz, nz);
            
            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.charge_density = (SpinResolved_DoubleMatrix)input_dict["SpinResolved_Density"]; else this.charge_density = new SpinResolved_DoubleMatrix(ny, nz);
            if (input_dict.ContainsKey("Potential")) this.band_offset = new Band_Data((DoubleMatrix)input_dict["Potential"]); else this.band_offset = new Band_Data(new DoubleMatrix(ny, nz));
            
            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            Console.WriteLine("Creating density calculation class");
            // create classes and initialise
            TwoD_ThomasFermiSolver dens_solv = new TwoD_ThomasFermiSolver(band_structure, acceptor_conc, donor_conc, acceptor_energy, donor_energy, temperature, dy, dz, ny, nz);
            double bottom_bc = dens_solv.Get_Chemical_Potential(nz - 1);
            Console.WriteLine("Density calculation class completed");

            Console.WriteLine("Creating potential solving class");
            TwoD_PoissonSolver pois_solv = new TwoD_PoissonSolver(dy, dz, ny, nz, bottom_bc, layer_depths, using_flexPDE, flexPDE_input, flexPDE_location, freeze_dopents, tol);
            Console.WriteLine("Potential solving class completed");

            // initialise dopents and band energies
            DoubleMatrix dopent_charge_density;
            if (freeze_dopents)
            {
                dopent_charge_density = Create_Dopent_Density_File(pois_solv).mat;
                // and use this density to calculate the initial band offset
                Console.WriteLine("Calculating initial band offset using dopent density");
                band_offset = new Band_Data(pois_solv.Get_Band_Energy(dopent_charge_density));
            }
            else
            {
                // but if there is no such density, start from no charge
                Console.WriteLine("Calculating initial band offset from scratch");
                band_offset = new Band_Data(pois_solv.Get_Band_Energy(new DoubleMatrix(ny, nz, 0.0)));
            }

            // run self consistent loop
            Console.WriteLine("Starting self-consistent density-potential solver loop");
            int count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString() + "\tConvergence factor:\t" + pois_solv.Convergence_Factor.ToString());

                // calculate the total charge density for this band offset
                charge_density = dens_solv.Get_TwoD_ChargeDensity(band_offset.mat);

                // solve the band energy for the given charge density and mix in with the old band energy
                Band_Data new_band_energy = new Band_Data(pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Matrix));
                pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                // change the potential mixing parameter
                //if ((count + 1) % potential_mixing_rate == 0)
                //    alpha = pois_solv.Renew_Mixing_Parameter(potential, new_potential, alpha_prime, alpha);

                // transfer new potential array to potential array
                //potential = new_potential;
                count++;
            }

            dens_solv.Output(charge_density, "charge_density.dat");
            pois_solv.Output(Input_Band_Structure.Expand_BandStructure(band_structure, ny) / 2.0 - band_offset, "potential.dat");
        }

        /// <summary>
        /// Finds whether there is a dopent density file created by the 1D band structure calculation and expands it into a 2D band_structure
        /// </summary>
        private Band_Data Create_Dopent_Density_File(TwoD_PoissonSolver pois_solv)
        {
            string[] raw_data;
            DoubleVector density;

            // check for dens.dat files
            if (File.Exists("dens_1D.dat"))
            {
                raw_data = File.ReadAllLines("dens_1D.dat");

                // separate the data into coordinate and density parts
                int no_points = int.Parse(raw_data[0].Split(' ')[1]);
                // take this data and save to density vector
                density = new DoubleVector(no_points);
                for (int i = 0; i < no_points; i++)
                    density[i] = float.Parse(raw_data[no_points + i + 3]);

                // check that the lattice spacing in the growth direction is the same
                double dz_data = double.Parse(raw_data[2]) - double.Parse(raw_data[1]);
                if (Math.Abs(dz - dz_data) > tol)
                    throw new Exception("Error - dz from input file is " + dz.ToString() + " whilst dz from 1d density is " + dz_data.ToString());
            }
            else
                throw new Exception("Error - trying to save out a 1D band structure to 2D failed... maybe the system is too cold?");

            // save the density into a specific file for the frozen dopent density
            Band_Data expanded_density = Input_Band_Structure.Expand_BandStructure(density, ny);
            pois_solv.Save_Density(expanded_density, "dens_2D.dat");
            pois_solv.Save_Density(expanded_density, "dopents_frozen_dens_2D.dat");

            return expanded_density;
        }
    }
}
