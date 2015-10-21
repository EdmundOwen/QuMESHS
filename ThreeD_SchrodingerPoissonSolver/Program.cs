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
            CenterSpace.NMath.Core.NMathConfiguration.LicenseKey = License.NMath_License_Key;

            Console.WriteLine("Program starting");

            Console.WriteLine("Loading input parameters from file");
            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Solver_Config.txt");
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

            // temporarily, just set the output suffix to something boring
            inputs.Add("output_suffix", ".dat");

            // initialise the band structure experiment
            Experiment exp = new Experiment();
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();

            Console.WriteLine("Performing density dopent calculation");
            Dictionary<string, object> inputs_init = new Dictionary<string, object>();
            inputs_init = inputs.Where(s => s.Key.ToLower().EndsWith("_1d")).ToDictionary(dict => dict.Key.Remove(dict.Key.Length - 3), dict => dict.Value);
            inputs_init.Add("BandStructure_File", inputs["BandStructure_File"]);
            inputs_init.Add("T", inputs["T"]);
            inputs_init.Add("output_suffix", "_1d.dat");

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
            tmp_dop_dens.Save_3D_Data("dens_3D_dopents.dat", (double)inputs["dx"] * ((double)inputs["nx"] + 1.0) / ((double)inputs["nx_1d"] - 1.0), y_scaling * (double)inputs["dy"] * ((double)inputs["ny"] + 1.0) / ((double)inputs["ny_1d"] - 1.0),  z_scaling * (double)inputs["dz_1d"], -1.0 * (double)inputs["dx"] * ((double)inputs["nx"] + 1.0) / 2.0, -1.0 * y_scaling * (double)inputs["dy"] * ((double)inputs["ny"] + 1.0) / 2.0, z_scaling * dopent_layer.Zmin);
            Console.WriteLine("Saved 1D dopent density");
            
            Console.WriteLine("Starting experiment");
            exp.Initialise(inputs);
            Console.WriteLine("Experiment initialised");
            exp.Run();
            Console.WriteLine("Experiment complete");
        }
    }
}
