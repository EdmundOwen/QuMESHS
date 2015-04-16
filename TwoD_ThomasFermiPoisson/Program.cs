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
            Console.WriteLine("Setting Centerspace key");

            // set nmath license key

            Console.WriteLine("Program starting");

            Console.WriteLine("Loading input parameters from file");
            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");
            Console.WriteLine("Input parameters loaded");

            ////////////////////////////////////////////////
            //
            //    EDITS FOR BATCH RUNS
            //
            ////////////////////////////////////////////////

            // read in the value of vsg to be used
            Console.WriteLine("Enter split gate voltage");
   //         inputs["split_V"] = double.Parse(Console.ReadLine());
            int index = int.Parse(args[0]);

            int maxval = 9;
            int i1 = index % maxval;
            int i2 = (index - i1) / maxval;

            // split gate with bias
//            double v1 = -0.5 - 0.25 * (double)i1;
//            double v2 = -0.5 - 0.02 * (double)i2;
//            if (v1 + v2 < -2.2 || v1 < v2)
//                return;
//            inputs["split_V1"] = v1;
//            inputs["split_V2"] = v2;
//            Console.WriteLine("Setting \"split_V1\" to " + ((double)inputs["split_V1"]).ToString() + "V");
//            Console.WriteLine("Setting \"split_V2\" to " + ((double)inputs["split_V2"]).ToString() + "V");
//            inputs["top_V"] = 0.0;
//            Console.WriteLine("Setting \"top_V\" to " + ((double)inputs["top_V"]).ToString() + "V");
//            inputs["output_suffix"] = "_sg1" + ((double)inputs["split_V1"]).ToString("F2") + "_sg2" + ((double)inputs["split_V2"]).ToString("F2") + ".dat";

            //top gated with constant side gate
            inputs["split_V"] = -0.5;// -0.46 - 0.005 * (double)i1;
            Console.WriteLine("Setting \"split_V\" to " + ((double)inputs["split_V"]).ToString() + "V");
            inputs["top_V"] = -1.65;// -1.7 + -0.01 * (double)i2;
            Console.WriteLine("Setting \"top_V\" to " + ((double)inputs["top_V"]).ToString() + "V");
            inputs["output_suffix"] = "_sg" + ((double)inputs["split_V"]).ToString("F3") + "_tg" + ((double)inputs["top_V"]).ToString("F3") + ".dat";

            ////////////////////////////////////////////////

            // check to make sure it's negative
            if ((double)inputs["split_V"] > 0)
            {
                Console.WriteLine("\"split_V\" has been set positive at " + ((double)inputs["split_V"]).ToString() + "V.  Are you sure you want to do this?");
                Console.ReadKey();
            }

            // initialise the band structure experiment
            Experiment exp = new Experiment();
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();

            // check if we should start from a precalculated density
            // consistency of band-structure, etc is the responsibility of the user...
            //if (!(bool)inputs["hot_start"])
            {
                Console.WriteLine("Performing density dopent calculation");
                Dictionary<string, object> inputs_init = new Dictionary<string, object>();
                Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
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

                //Input_Band_Structure.Expand_BandStructure(exp_init.Dopent_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D_dopents.dat", (double)inputs["dy"] * (double)inputs["ny"] / (double)inputs_init["ny_1d"], (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * (double)inputs["ny"] / 2.0, Geom_Tool.Get_Zmin(exp_init.Layers));
               
                // this is a scaled version for the dopents!
                double scaling_factor = ((double)inputs["ny"] * (double)inputs["dy"]) / ((double)inputs["nz"] * (double)inputs["dz"]);
                Input_Band_Structure.Expand_BandStructure(exp_init.Dopent_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D_dopents.dat", (double)inputs["dy"] * ((double)inputs["ny"] + 2.0) / ((double)inputs_init["ny_1d"] - 1.0), scaling_factor * (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * ((double)inputs["ny"] + 2.0) / 2.0, scaling_factor * Geom_Tool.Get_Zmin(exp_init.Layers));
         //       Input_Band_Structure.Expand_BandStructure(exp_init.Carrier_Density, (int)(double)inputs_init["ny_1d"]).Spin_Summed_Data.Save_2D_Data("dens_2D.dat", (double)inputs["dy"] * ((double)inputs["ny"] + 2.0) / ((double)inputs_init["ny_1d"] - 1.0), (double)inputs_init["dz"], -1.0 * (double)inputs["dy"] * ((double)inputs["ny"] + 2.0) / 2.0, Geom_Tool.Get_Zmin(exp_init.Layers));
                Console.WriteLine("Saved 1D dopent density");
            }

            if ((bool)inputs["batch_run"])
            {
                Run_Multiple_SGs(inputs);
            }
            else
            {
                Console.WriteLine("Starting experiment");
                exp.Initialise_Experiment(inputs);
                // check that the dz_pot are the same for both simulations as this is needed for the interpolation of SpinResolved_Density
                if (!(bool)inputs["hot_start"] && exp_init.Dz_Pot != exp.Dz_Pot)
                    throw new Exception("Error - the dz values for the potentials must be the same for \"Input_Parameters.txt\" and \"Input_Parameters_1D.txt\"");
                Console.WriteLine("Experiment initialised");
                exp.Run();

                Console.WriteLine("Experiment complete");
            }
        }

        static void Rename_Results(double split_V, double top_V)
        {
            Console.WriteLine();
            Console.WriteLine("Renaming files");
            File.Copy("bare_pot.dat", "bare_pot_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            File.Copy("potential.dat", "pot_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            File.Copy("dens_2D_up_raw.dat", "dens_2D_up_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            File.Copy("dens_2D_down_raw.dat", "dens_2D_down_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            File.Copy("energies.dat", "energies_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            File.Copy("xc_pot.dat", "xc_pot_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            File.Copy("pot_KS.dat", "pot_KS_sg" + split_V.ToString("F2") + "_tg" + top_V.ToString("F2") + ".dat", true);
            Console.WriteLine("Files renamed");
        }

        static void Run_Multiple_TGs(TwoD_ThomasFermiPoisson.Experiment exp, Dictionary<string, object> dict)
        {
            int init = (int)Math.Abs((int)(double)dict["init_val"]);
            int final = (int)Math.Abs((int)(double)dict["final_val"]);
            double interval = (double)dict["interval"] * Math.Sign((double)dict["final_val"] - (double)dict["init_val"]);

            for (int i = init; i < final; i++)
            {
                double tg = i * interval;
                dict["top_V"] = tg;
                exp.Initialise_Experiment(dict);
                Console.WriteLine("Experiment initialised for tg = " + tg.ToString() + "V");
                exp.Run();
                File.Copy("dens_2D_up_raw.dat", "dens_2D_up_sg05_tg" + i.ToString("000") + ".dat");
                File.Copy("dens_2D_down_raw.dat", "dens_2D_down_sg05_tg" + i.ToString("000") + ".dat");
                File.Copy("energies.dat", "energies_sg05_tg" + i.ToString("000") + ".dat");
                Console.WriteLine("Experiment complete");
            }
        }

        static void Anneal_Multiple_TGs(TwoD_ThomasFermiPoisson.Experiment exp, Dictionary<string, object> dict)
        {
            int init = (int)Math.Abs((int)(double)dict["init_val"]);
            int final = (int)Math.Abs((int)(double)dict["final_val"]);
            double interval = (double)dict["interval"] * Math.Sign((double)dict["final_val"] - (double)dict["init_val"]);

            for (int i = init; i < final; i++)
            {
                File.Copy("..//dens_2D_up_sg05_tg" + i.ToString("000") + ".dat", "dens_2d_up_raw.dat", true);
                File.Copy("..//dens_2D_down_sg05_tg" + i.ToString("000") + ".dat", "dens_2d_down_raw.dat", true);

                double tg = i * interval;
                dict["top_V"] = tg;
                exp.Initialise_Experiment(dict);
                Console.WriteLine("Experiment initialised for tg = " + tg.ToString() + "V");
                exp.Run();
                File.Copy("dens_2D_up_raw.dat", "dens_2D_up_sg05_tg" + i.ToString("000") + ".dat");
                File.Copy("dens_2D_down_raw.dat", "dens_2D_down_sg05_tg" + i.ToString("000") + ".dat");
                File.Copy("energies.dat", "energies_sg05_tg" + i.ToString("000") + ".dat");
                Console.WriteLine("Experiment complete");
            }
        }

        static void Run_Multiple_SGs(Dictionary<string, object> dict)
        {
            TwoD_ThomasFermiPoisson.Experiment exp;

            int init = (int)(double)dict["init_val"];
            int final = (int)(double)dict["final_val"];
            double interval = (double)dict["interval"];

            // and set an integer value for the split gate voltage.... this will be used in the save file names
            double split_V_init = (double)dict["split_V"];

            // set the dictionary values equal to the values we're going to use in the simulation loop
//            int vtg_start = init;
//            int vtg_end = final;
                        
            int vtg_start = init;
            int vtg_end = final;
            int no_vsg = (int)(double)dict["no_Vsg"];
            double delta_vsg = (double)dict["delta_Vsg"];

            for (int i = 0; i < no_vsg; i++)
            {
                double split_V = split_V_init + (double)i * delta_vsg;

                // generate null raw data
                string[] tmp_data = new string[(int)(double)dict["ny_dens"] * (int)(double)dict["nz_dens"]];
                for (int k = 0; k < tmp_data.Length; k++)
                    tmp_data[k] = "0";
                // set starting file name
                dict["spin_up_file"] = "dens_2D_up_sg" + split_V.ToString("F3") + "_tg" + ((double)(vtg_start - 1) * interval).ToString("F3") + ".dat";
                dict["spin_down_file"] = "dens_2D_down_sg" + split_V.ToString("F3") + "_tg" + ((double)(vtg_start - 1) * interval).ToString("F3") + ".dat";
                // and write data there
                if (!File.Exists((string)dict["spin_up_file"]))
                    File.WriteAllLines((string)dict["spin_up_file"], tmp_data);
                if (!File.Exists((string)dict["spin_down_file"]))
                    File.WriteAllLines((string)dict["spin_down_file"], tmp_data);
                File.WriteAllLines((string)dict["surface_charge_file"], new string[] { ((double)dict["surface_charge"]).ToString() });

                for (int j = vtg_start; j < vtg_end + 1; j++)
                {
                    exp = new Experiment();
                    double top_V = (double)j * interval;
                    dict["top_V"] = top_V;
                    dict["split_V"] = split_V;
                    dict["spin_up_file"] = "dens_2D_up_sg" + split_V.ToString("F3") + "_tg" + ((double)(j - 1) * interval).ToString("F3") + ".dat";
                    dict["spin_down_file"] = "dens_2D_down_sg" + split_V.ToString("F3") + "_tg" + ((double)(j - 1) * interval).ToString("F3") + ".dat";

                    exp.Initialise_Experiment(dict);
                    Console.WriteLine("Experiment initialised for sg = " + ((double)dict["split_V"]).ToString() + "V, tg = " + ((double)dict["top_V"]).ToString() + "V");
                    exp.Run();

                    //                File.Copy("dens_2D_up_raw.dat", "dens_2D_up_sg" + i.ToString("0000") + "_tg" + ((double)dict["top_V"] * 100).ToString("0000") + ".dat", true);
                    //                File.Copy("dens_2D_down_raw.dat", "dens_2D_down_sg" + i.ToString("0000") + "_tg" + ((double)dict["top_V"] * 100).ToString("0000") + ".dat", true);
                    //                File.Copy("energies.dat", "energies_sg" + i.ToString("0000") + "_tg" + ((double)dict["top_V"] * 100).ToString("0000") + ".dat", true);
                    //                File.Copy("split_gate_final.pg6", "split_gate_final_sg" + i.ToString("0000") + "_tg" + ((double)dict["top_V"] * 100).ToString("0000") + ".pg6", true);

                    File.Copy("bare_pot.dat", "bare_pot_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("potential.dat", "pot_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("dens_2D_up_raw.dat", "dens_2D_up_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("dens_2D_down_raw.dat", "dens_2D_down_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("energies.dat", "energies_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("xc_pot.dat", "xc_pot_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("pot_KS.dat", "pot_KS_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("ks_ke.dat", "ks_ke_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".dat", true);
                    File.Copy("split_gate_final.pg6", "split_gate_final_sg" + split_V.ToString("F3") + "_tg" + top_V.ToString("F3") + ".pg6", true);
                    Console.WriteLine("Experiment complete");
                }
            }
        }
    }
}
