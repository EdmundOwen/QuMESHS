using System;
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

            Console.WriteLine("Performing density dopent calculation");
            Dictionary<string, object> inputs_init = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();
            exp_init.Initialise(inputs_init);
            exp_init.Run();
            inputs.Add("SpinResolved_Density", exp_init.Charge_Density);
            inputs.Add("Band_Offset", exp_init.Band_Offset);
            Console.WriteLine("Calculated 1D density for dopents");

            Input_Band_Structure.Expand_BandStructure(exp_init.Charge_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D_donors.dat", (double)inputs["dy"] * (double)inputs["ny"] / (double)inputs_init["ny_1d"], (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * (double)inputs["ny"] / 2.0, Geom_Tool.Get_Zmin(exp_init.Layers));
            Input_Band_Structure.Expand_BandStructure(exp_init.Charge_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D.dat", (double)inputs["dy"] * (double)inputs["ny"] / (double)inputs_init["ny_1d"], (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * (double)inputs["ny"] / 2.0, Geom_Tool.Get_Zmin(exp_init.Layers));
            Console.WriteLine("Saved 1D dopent density");

            Band_Data band_offset = exp_init.Band_Offset;
            ILayer[] layers = exp_init.Layers;
            OneD_ThomasFermiPoisson.OneD_PoissonSolver tmp_pois_solv = new OneD_ThomasFermiPoisson.OneD_PoissonSolver(exp_init, false, "", "", 0.0);
            inputs.Add("surface_charge", tmp_pois_solv.Get_Surface_Charge(band_offset, layers) * (double)inputs_init["dz"]);

            Console.WriteLine("Starting experiment");
            Experiment exp = new Experiment();
            Console.WriteLine("Created experiment");
            exp.Initialise_Experiment(inputs);
            exp.Run();
            Console.WriteLine("Experiment complete");
        }
    }
}
