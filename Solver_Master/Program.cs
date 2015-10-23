/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;
using Iterative_Greens_Function_Test;

namespace Solver_Master
{
    class Program
    {
        static void Main(string[] args)
        {
            // set nmath license key
            CenterSpace.NMath.Core.NMathConfiguration.LicenseKey = License.NMath_License_Key;

            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Solver_Config.txt");
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, args[1]);

            int no_runs = 1;
            bool batch_run = false;
            if (inputs.ContainsKey("batch_run"))
                batch_run = (bool)inputs["batch_run"];
            if (inputs.ContainsKey("no_runs"))
                no_runs = (int)(double)inputs["no_runs"];

            int first_run = 0;
            if (args.Length != 0)
                first_run = int.Parse(args[0]);

            if (!inputs.ContainsKey("dim"))
                throw new KeyNotFoundException("Error - you must define the dimensionality of the system you want to solve!");

            int dim = (int)(double)inputs["dim"];
            Calculate_1D_Band_Structure(inputs);

            for (int batch_no = first_run; batch_no < first_run + no_runs; batch_no++)
            {
                inputs["batch_no"] = batch_no;
                IExperiment exp;
                if (batch_run)
                    Prepare_Batch_Inputs(inputs, batch_no);

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

        static void Prepare_Batch_Inputs(Dictionary<string, object> inputs, int batch_no)
        {
            // remove any input parameters from previous runs
            inputs.Remove("output_suffix");

            // get the boundary conditions which will be editted in the batch
            string[] batch_params = ((string)inputs["batch_params"]).TrimStart('{').TrimEnd('}').Split(',');
            string output_suffix = "";

            int[] nos = new int[batch_params.Length];
            double[] deltas = new double[batch_params.Length];
            double[] inits = new double[batch_params.Length];
            int[] batch_index = new int[batch_params.Length];

            for (int i = 0; i < batch_params.Length; i++)
            {
                if (batch_params[i].Contains('='))
                    break;

                string tmp_bc = batch_params[i].Trim();

                // check that the batch bc is correctly labelled
                if (!inputs.ContainsKey(tmp_bc))
                    throw new KeyNotFoundException("Error - there should be a dummy option for the requested batch variable \"" + tmp_bc + "\"");
                if (!inputs.ContainsKey("no_" + tmp_bc))
                    throw new KeyNotFoundException("Error - there is no \"no_" + tmp_bc + "\" option defined in the input file");
                if (!inputs.ContainsKey("delta_" + tmp_bc))
                    throw new KeyNotFoundException("Error - there is no \"delta_" + tmp_bc + "\" option defined in the input file");
                if (!inputs.ContainsKey("init_" + tmp_bc))
                    throw new KeyNotFoundException("Error - there is no \"init_" + tmp_bc + "\" option defined in the input file");

                // get the batch information
                nos[i] = (int)(double)inputs["no_" + tmp_bc];
                deltas[i] = (double)inputs["delta_" + tmp_bc];
                inits[i] = (double)inputs["init_" + tmp_bc];
            }

            // process the batch number
            batch_index[0] = batch_no % nos[0];
            for (int i = 1; i < batch_params.Length; i++)
            {
                int mod_val = nos[0];
                for (int j = 1; j < i; j++)
                    mod_val *= nos[j];

                batch_index[i] = (batch_no - batch_no % mod_val) / mod_val;
            }

            // change the input library according to the requested batch input
            for (int i = 0; i < batch_params.Length; i++)
                if (!batch_params[i].Contains('='))
                {
                    inputs[batch_params[i].Trim()] = inits[i] + deltas[i] * (double)batch_index[i];
                    Console.WriteLine("Batch parameter \"" + batch_params[i].Trim() + "\" set to " + ((double)inputs[batch_params[i].Trim()]).ToString());
                }
                else
                {
                    string[] vals = batch_params[i].Split('=');
                    inputs[vals[0].Trim()] = (double)inputs[vals[1].Trim()];
                }

            // reset the "voltages" value to reflect the new batch values
            int count = 0;
            string voltages = "{";
            while (inputs.ContainsKey("V" + count.ToString()))
            {
                voltages += ((double)inputs["V" + count.ToString()]).ToString() + ", ";
                count++;
            }
            inputs["voltages"] = voltages.Remove(voltages.Length - 2) + "}";

            // and generate the relevant batch suffix
            for (int i = 0; i < batch_params.Length; i++)
                if (!batch_params[i].Contains('='))
                    output_suffix += "_" + batch_params[i].Trim() + "_" + ((double)inputs[batch_params[i].Trim()]).ToString("F3");
            inputs.Add("output_suffix", output_suffix + ".dat");
        }

        static void Calculate_1D_Band_Structure(Dictionary<string, object> inputs)
        {
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();

            Console.WriteLine("Performing density dopent calculation");
            Dictionary<string, object> inputs_init = new Dictionary<string, object>();
            if ((int)(double)inputs["dim"] != 1)
            {
                inputs_init = inputs.Where(s => s.Key.ToLower().EndsWith("_1d")).ToDictionary(dict => dict.Key.Remove(dict.Key.Length - 3), dict => dict.Value);
                inputs_init.Add("BandStructure_File", inputs["BandStructure_File"]);
                inputs_init.Add("output_suffix", "_1d.dat");
                inputs_init.Add("T", inputs["T"]);
            }
            else
                inputs_init = inputs.Where(s => s.Key.ToLower().EndsWith("")).ToDictionary(dict => dict.Key, dict => dict.Value);


            //    Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
            exp_init.Initialise(inputs_init);
            exp_init.Run();
            inputs.Add("Carrier_Density", exp_init.Carrier_Density);
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
