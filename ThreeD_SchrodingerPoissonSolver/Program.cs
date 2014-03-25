using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;

namespace ThreeD_SchrodingerPoissonSolver
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
            inputs.Add("SpinResolved_Density", exp_init.Carrier_Density);
            inputs.Add("Band_Offset", exp_init.Chemical_Potential);
            Console.WriteLine("Calculated 1D density for dopents");

            Band_Data chem_pot = exp_init.Chemical_Potential;
            ILayer[] layers = exp_init.Layers;
            OneD_ThomasFermiPoisson.OneD_PoissonSolver tmp_pois_solv = new OneD_ThomasFermiPoisson.OneD_PoissonSolver(exp_init, false, "", "", 0.0);
            inputs.Add("surface_charge", tmp_pois_solv.Get_Surface_Charge(chem_pot, layers) * (double)inputs_init["dz"]);
            inputs.Add("bottom_bc", exp_init.Bottom_BC);
            
            Console.WriteLine("Starting experiment");
            Experiment exp = new Experiment();
            Console.WriteLine("Created experiment");
            exp.Initialise_Experiment(inputs);

            int nx_init = (int)(double)inputs_init["nx_1d"]; int ny_init = (int)(double)inputs_init["ny_1d"];
            double dx_init = exp.Nx_Pot * exp.Dx_Pot / nx_init; double dy_init = exp.Ny_Pot * exp.Dy_Pot / ny_init;
            
            Input_Band_Structure.Expand_BandStructure(exp_init.Dopent_Density, nx_init, ny_init).Spin_Summed_Data.Save_3D_Data("dens_3D_dopents.dat", dx_init, dy_init, exp.Dz_Pot, exp.Xmin_Pot, exp.Ymin_Pot, exp.Zmin_Pot);
            Input_Band_Structure.Expand_BandStructure(exp_init.Carrier_Density, nx_init, ny_init).Spin_Summed_Data.Save_3D_Data("dens_3D.dat", dx_init, dy_init, exp.Dz_Pot, exp.Xmin_Pot, exp.Ymin_Pot, exp.Zmin_Pot);
            Console.WriteLine("Saved 1D dopent density");

            exp.Run();
            Console.WriteLine("Experiment complete");
        }
    }
}
