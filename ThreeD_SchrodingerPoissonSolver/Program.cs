using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace ThreeD_SchrodingerPoissonSolver
{
    class Program
    {
        static void Main(string[] args)
        {
            // set nmath license key

            Console.BufferHeight = Int16.MaxValue - 1;
            Console.WriteLine("Program starting");

            Console.WriteLine("Loading input parameters from file");
            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");
            Console.WriteLine("Input parameters loaded");

            // read in the value of vsg to be used
            Console.WriteLine("Enter split gate voltage");
            inputs["split_V"] = double.Parse(Console.ReadLine());
            Console.WriteLine("Setting \"split_V\" to " + ((double)inputs["split_V"]).ToString() + "V");

            // check to make sure it's negative
            if ((double)inputs["split_V"] > 0)
            {
                Console.WriteLine("\"split_V\" has been set positive at " + ((double)inputs["split_V"]).ToString() + "V.  Are you sure you want to do this?");
                Console.ReadKey();
            }

            // initialise the band structure experiment
            Experiment exp = new Experiment();
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();

            Console.WriteLine("Performing density dopent calculation");
            Dictionary<string, object> inputs_init = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
            exp_init.Initialise(inputs_init);
            exp_init.Run();
            inputs.Add("SpinResolved_Density", exp_init.Carrier_Density);
            inputs.Add("Chemical_Potential", exp_init.Chemical_Potential);
            // get the frozen out surface charge at 70K
            if (!inputs.ContainsKey("surface_charge")) inputs.Add("surface_charge", exp_init.Surface_Charge(70.0));
            else Console.WriteLine("Surface charge set from Input_Parameters.txt to " + ((double)inputs["surface_charge"]).ToString());
            Console.WriteLine("Calculated 1D density for dopents");

            int nx_init = (int)(double)inputs_init["nx_1d"]; int ny_init = (int)(double)inputs_init["ny_1d"];
            double dx_init = (double)inputs["nx"] * (double)inputs["dx"] / nx_init; double dy_init = (double)inputs["ny"] * (double)inputs["ny"] / ny_init;

            // this is a scaled version for the dopents!
            double y_scaling = ((double)inputs["nx"] * (double)inputs["dx"]) / ((double)inputs["ny"] * (double)inputs["dy"]);
            double z_scaling = ((double)inputs["nx"] * (double)inputs["dx"]) / ((double)inputs["nz"] * (double)inputs["dz"]);
            Input_Band_Structure.Expand_BandStructure(exp_init.Dopent_Density, nx_init, ny_init).Spin_Summed_Data.Save_3D_Data("dens_3D_dopents.dat", dx_init, dy_init * y_scaling, (double)inputs["dz"] * z_scaling, (double)inputs["xmin_pot"], (double)inputs["ymin_pot"] * y_scaling, Geom_Tool.Get_Zmin(exp_init.Layers) * z_scaling);
            Input_Band_Structure.Expand_BandStructure(exp_init.Carrier_Density, nx_init, ny_init).Spin_Summed_Data.Save_3D_Data("dens_3D.dat", dx_init, dy_init * y_scaling, (double)inputs["dz"] * z_scaling, (double)inputs["xmin_pot"], (double)inputs["ymin_pot"] * y_scaling, Geom_Tool.Get_Zmin(exp_init.Layers) * z_scaling);
            Console.WriteLine("Saved 1D dopent density");
            
            Console.WriteLine("Starting experiment");
            exp.Initialise_Experiment(inputs);
            Console.WriteLine("Experiment initialised");
            exp.Run();
            Console.WriteLine("Experiment complete");
        }
    }
}
