using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;

namespace Solver_Master
{
    class Program
    {
        static void Main(string[] args)
        {
            // set nmath license key

            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Solver_Config.txt");
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");

            int no_runs = 1;
            if (inputs.ContainsKey("batch_run"))
                no_runs = (int)(double)inputs["no_runs"];

            int dim = (int)(double)inputs["dim"];
            if (dim != 1)
                Calculate_1D_Band_Structure(inputs);

            for (int i = 0; i < no_runs; i++)
            {
                IExperiment exp;
                if (no_runs != 1)
                    Prepare_Batch_Inputs(inputs, i);

                switch (dim)
                {
                    case 1:
                        exp = new OneD_ThomasFermiPoisson.Experiment();
                        break;
                    case 2:
                        exp = new TwoD_ThomasFermiPoisson.Experiment();
                        break;
                    case 3:
                        exp = new ThreeD_SchrodingerPoissonSolver.Experiment();
                        break;
                    default:
                        throw new ArgumentException("Error - Requested dimension for the system dim = " + dim.ToString() + " not supported!");
                }

                exp.Initialise(inputs);
                exp.Run();
            }
        }

        private static void Prepare_Batch_Inputs(Dictionary<string, object> inputs, int i)
        {
            throw new NotImplementedException();
        }

        static void Calculate_1D_Band_Structure(Dictionary<string, object> inputs)
        {
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();

            Console.WriteLine("Performing density dopent calculation");
            Dictionary<string, object> inputs_init = new Dictionary<string, object>();
            inputs_init = inputs.Where(s => s.Key.ToLower().EndsWith("_1d")).ToDictionary(dict => dict.Key.Remove(dict.Key.Length - 3), dict => dict.Value);
            inputs_init.Add("BandStructure_File", inputs["BandStructure_File"]);
            inputs_init.Add("T", inputs["T"]);
            inputs_init.Add("output_suffix", "_1d.dat");

            //    Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
            exp_init.Initialise(inputs_init);
            exp_init.Run();
            inputs.Add("SpinResolved_Density", exp_init.Carrier_Density);
            inputs.Add("Dopent_Density", exp_init.Dopent_Density);
            inputs.Add("Chemical_Potential", exp_init.Chemical_Potential);
            inputs.Add("nz_pot_1d", inputs_init["nz"]);
            inputs.Add("zmin_pot_1d", inputs_init["zmin"]);
            inputs.Add("zmax_pot_1d", inputs_init["zmax"]);
            // get the frozen out surface charge at 70K
            if (!inputs.ContainsKey("surface_charge")) inputs.Add("surface_charge", exp_init.Surface_Charge(70.0));
            else Console.WriteLine("Surface charge set from Input_Parameters.txt to " + ((double)inputs["surface_charge"]).ToString());
            Console.WriteLine("Calculated 1D density for dopents");

            if ((int)(double)inputs["dim"] == 2)
            {
                // create a scaled data file containing the dopent density
                double scaling_factor = ((double)inputs["ny"] * (double)inputs["dy"]) / ((double)inputs["nz"] * (double)inputs["dz"]);
                Input_Band_Structure.Expand_BandStructure(exp_init.Dopent_Density, (int)(double)inputs["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D_dopents.dat", (double)inputs["dy"] * ((double)inputs["ny"] + 2.0) / ((double)inputs["ny_1d"] - 1.0), scaling_factor * (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * ((double)inputs["ny"] + 2.0) / 2.0, scaling_factor * Geom_Tool.Get_Zmin(exp_init.Layers));
            }
            else if ((int)(double)inputs["dim"] == 3)
            {
                // this is a scaled version for the dopents!
                double y_scaling = ((double)inputs["nx"] * (double)inputs["dx"]) / ((double)inputs["ny"] * (double)inputs["dy"]);
                double z_scaling = ((double)inputs["nx"] * (double)inputs["dx"]) / ((double)inputs["nz"] * (double)inputs["dz"]);
                // extract the dopent layer (leaving the top and bottom set to zero)
                int dopent_min = -1;
                int dopent_max = -2;
                ILayer dopent_layer = exp_init.Layers[0];
                for (int i = 0; i < exp_init.Layers.Length; i++)
                    if (exp_init.Layers[i].Donor_Conc != 0.0 || exp_init.Layers[i].Acceptor_Conc != 0)
                    {
                        dopent_layer = exp_init.Layers[i];
                        dopent_min = (int)Math.Round((dopent_layer.Zmin - Geom_Tool.Get_Zmin(exp_init.Layers)) / (int)(double)inputs_init["dz"]);
                        dopent_max = (int)Math.Round((dopent_layer.Zmax - Geom_Tool.Get_Zmin(exp_init.Layers)) / (int)(double)inputs_init["dz"]);
                    }
                Band_Data tmp_dop_dens_1D = new Band_Data(dopent_max - dopent_min, 0.0);
                for (int i = dopent_min + 1; i < dopent_max - 1; i++)
                    tmp_dop_dens_1D.vec[i - dopent_min] = exp_init.Dopent_Density.Spin_Summed_Data.vec[i];
                // and expand into the correct data structure
                Band_Data tmp_dop_dens = Input_Band_Structure.Expand_BandStructure(tmp_dop_dens_1D.vec, (int)(double)inputs["nx_1d"], (int)(double)inputs["ny_1d"]);
                tmp_dop_dens.Save_3D_Data("dens_3D_dopents.dat", (double)inputs["dx"] * ((double)inputs["nx"] + 1.0) / ((double)inputs["nx_1d"] - 1.0), y_scaling * (double)inputs["dy"] * ((double)inputs["ny"] + 1.0) / ((double)inputs["ny_1d"] - 1.0), z_scaling * (double)inputs["dz_1d"], -1.0 * (double)inputs["dx"] * ((double)inputs["nx"] + 1.0) / 2.0, -1.0 * y_scaling * (double)inputs["dy"] * ((double)inputs["ny"] + 1.0) / 2.0, z_scaling * dopent_layer.Zmin);
                Console.WriteLine("Saved 1D dopent density");
            }
        }
    }
}
