﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace TwoD_ThomasFermiPoisson
{
    class Program
    {
        static void Main(string[] args)
        {
            // set nmath license key

            Console.WriteLine("Program starting");

            Console.WriteLine("Loading input parameters from file");
            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");
            Console.WriteLine("Input parameters loaded");

            Experiment exp = new Experiment();
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();

            // check if we should start from a precalculated density
            // consistency of band-structure, etc is the responsibility of the user...
            if (!bool.Parse((string)inputs["hot_start"]))
            {
                Console.WriteLine("Performing density dopent calculation");
                Dictionary<string, object> inputs_init = new Dictionary<string, object>();
                Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
                exp_init.Initialise(inputs_init);
                exp_init.Run();
                inputs.Add("SpinResolved_Density", exp_init.Carrier_Density);
                inputs.Add("Chemical_Potential", exp_init.Chemical_Potential);
                Console.WriteLine("Calculated 1D density for dopents");

                Input_Band_Structure.Expand_BandStructure(exp_init.Dopent_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D_dopents.dat", (double)inputs["dy"] * (double)inputs["ny"] / (double)inputs_init["ny_1d"], (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * (double)inputs["ny"] / 2.0, Geom_Tool.Get_Zmin(exp_init.Layers));
                Input_Band_Structure.Expand_BandStructure(exp_init.Carrier_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D.dat", (double)inputs["dy"] * (double)inputs["ny"] / (double)inputs_init["ny_1d"], (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * (double)inputs["ny"] / 2.0, Geom_Tool.Get_Zmin(exp_init.Layers));
                Console.WriteLine("Saved 1D dopent density");

                Band_Data band_offset = exp_init.Chemical_Potential;
                ILayer[] layers = exp_init.Layers;
                OneD_ThomasFermiPoisson.OneD_PoissonSolver tmp_pois_solv = new OneD_ThomasFermiPoisson.OneD_PoissonSolver(exp_init, false, "", "", 0.0);
                inputs.Add("surface_charge", tmp_pois_solv.Get_Surface_Charge(band_offset, layers) * (double)inputs_init["dz"]);
            }

            
           // Console.WriteLine("Starting experiment");
           // exp.Initialise_Experiment(inputs);
           // // check that the dz_pot are the same for both simulations as this is needed for the interpolation of SpinResolved_Density
           // if (!bool.Parse((string)inputs["hot_start"]) && exp_init.Dz_Pot != exp.Dz_Pot)
           //     throw new Exception("Error - the dz values for the potentials must be the same for \"Input_Parameters.txt\" and \"Input_Parameters_1D.txt\"");
           // Console.WriteLine("Experiment initialised");
           // exp.Run();
           // Console.WriteLine("Experiment complete");

            Run_Multiple_TGs(exp, inputs);
        }

        static void Run_Multiple_TGs(TwoD_ThomasFermiPoisson.Experiment exp, Dictionary<string, object> dict)
        {
            for (int i = 1; i < 400; i++)
            {
                double tg = i * -0.01;
                dict["top_V"] = tg;
                exp.Initialise_Experiment(dict);
                Console.WriteLine("Experiment initialised for tg = " + tg.ToString() + "V");
                exp.Run();
                File.Copy("dens_2D_up_raw.dat", "dens_2D_up_sg13_tg" + i.ToString("000") + ".dat");
                File.Copy("dens_2D_down_raw.dat", "dens_2D_down_sg13_tg" + i.ToString("000") + ".dat");
                File.Copy("energies.dat", "energies_sg13_tg" + i.ToString("000") + ".dat");
                Console.WriteLine("Experiment complete");
            }
        }

        static void Run_Multiple_SGs(TwoD_ThomasFermiPoisson.Experiment exp, Dictionary<string, object> dict)
        {
            for (int i = 131; i < 400; i++)
            {
                double sg = i * -0.01;
                dict["split_V"] = sg;
                dict["spin_up_file"] = "dens_2D_up_sg" + i.ToString("000") + "_tgnat.dat";
                dict["spin_down_file"] = "dens_2D_down_sg" + i.ToString("000") + "_tgnat.dat";
                
                exp.Initialise_Experiment(dict);
                Console.WriteLine("Experiment initialised for sg = " + sg.ToString() + "V");
                exp.Run();
                //File.Copy("dens_2D_up.dat", "dens_2D_up_sg" + i.ToString("000") + "_tgnat.dat");
                //File.Copy("dens_2D_down.dat", "dens_2D_down_sg" + i.ToString("000") + "_tgnat.dat");
                File.Copy("energies.dat", "energies_sg" + i.ToString("000") + "_tgnat.dat");
                Console.WriteLine("Experiment complete");
            }
        }
    }
}
