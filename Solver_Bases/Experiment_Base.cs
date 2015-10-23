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
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace Solver_Bases
{
    public abstract class Experiment_Base : IExperiment
    {
        protected Dictionary<string, double> device_dimensions, boundary_conditions;

        protected SpinResolved_Data carrier_charge_density;
        protected SpinResolved_Data dopent_charge_density;
        protected SpinResolved_Data carrier_charge_density_deriv;
        protected SpinResolved_Data dopent_charge_density_deriv;
        protected Band_Data chem_pot;
        protected Band_Data x;
        protected Band_Data gphi;
        protected ILayer[] layers;

        // damping parameter for Bank-Rose convergence method
        protected double t = 0.5;

        // parameters for the density domain
        protected double dx_dens = 1.0, dy_dens = 1.0, dz_dens = 1.0;
        protected double xmin_dens = -1.0, ymin_dens = -1.0, zmin_dens = -1.0;
        protected int nx_dens = 1, ny_dens = 1, nz_dens = 1;

        // parameters for the potential domain
        protected double dx_pot, dy_pot, dz_pot;
        protected double xmin_pot = -1.0, ymin_pot = -1.0, zmin_pot = -1.0;
        protected int nx_pot = 1, ny_pot = 1, nz_pot = 1;

        protected double tol, tol_anneal;

        protected bool using_flexPDE = false;
        protected bool using_dealii = false;

        protected double initial_temperature = 300.0;       // initial temperature
        protected double current_temperature = 300.0;       // temperature currently being simulated
        protected double temperature;                       // target temperature

        protected int max_iterations = 1000;       // maximum number of runs before giving up and just outputting the data with a "not_converged" file flag
        protected bool converged = false;

        protected bool initial_run = true;      // should we do an initial run which stops before convergence in order to get a better initial density distribution
        protected int initial_run_steps = 5;    // and if yes, how many steps should it iterate for?

        protected bool no_dft = false;       // do not run with dft potential (ie. Hartree approximation)
        protected double dft_mixing_parameter = 0.3;    // Default mixing parameter for DFT
        protected bool hot_start = false;     // am the program starting from a precalculated density or do i start from scratch
        protected bool initialise_from_restart = false;     //is the program starting from a restart? (normally false)
        protected bool initialise_with_1D_data = false;     // should the program use the data used to calculate the dopents in order to find the initial density for the higher dimensional structure?

        // suffix for data output to identify between different runs... initially empty
        protected string output_suffix = "";

        public virtual void Initialise(Dictionary<string, object> input_dict)
        {
            // solver inputs
            Get_From_Dictionary<double>(input_dict, "tolerance", ref tol);
            tol_anneal = tol;
            Get_From_Dictionary<double>(input_dict, "anneal_tolerance", ref tol_anneal, true);
            Get_From_Dictionary(input_dict, "max_iterations", ref max_iterations, true);

            Get_From_Dictionary<bool>(input_dict, "initial_run", ref initial_run, true);
            Get_From_Dictionary(input_dict, "initial_run_steps", ref initial_run_steps, true);

            // will not use FlexPDE unless told to
            if (input_dict.ContainsKey("use_FlexPDE")) this.using_flexPDE = (bool)input_dict["use_FlexPDE"]; else using_flexPDE = false;
            // and the same for dealii
            if (input_dict.ContainsKey("use_deal.II")) this.using_dealii = (bool)input_dict["use_deal.II"]; else using_dealii = false;

            // check to make sure that we haven't asked to use more than one external program to calculate the potential
            if (using_flexPDE == true && using_dealii == true)
                throw new Exception("Error - Cannot use both FlexPDE and deal.II to calculate the potential!");

            // physical inputs
            Get_From_Dictionary<double>(input_dict, "init_T", ref initial_temperature, true);
            Get_From_Dictionary<double>(input_dict, "T", ref temperature);

            // get the band structure
            if (input_dict.ContainsKey("Layers"))
                layers = (ILayer[])input_dict["Layers"];
            else
            {
                if (input_dict.ContainsKey("BandStructure_File"))
                {
                    layers = Input_Band_Structure.Get_Layers((string)input_dict["BandStructure_File"]);
                    input_dict.Add("Layers", layers);
                }
                else throw new KeyNotFoundException("No band structure file found in input dictionary!");
            }

            // but try to get the specific values
            Get_From_Dictionary<double>(input_dict, "dx_dens", ref dx_dens, true);
            Get_From_Dictionary<double>(input_dict, "dy_dens", ref dy_dens, true);
            Get_From_Dictionary<double>(input_dict, "dz_dens", ref dz_dens, true);
            Get_From_Dictionary<double>(input_dict, "dx_pot", ref dx_pot, true);
            Get_From_Dictionary<double>(input_dict, "dy_pot", ref dy_pot, true);
            Get_From_Dictionary<double>(input_dict, "dz_pot", ref dz_pot, true);

            Get_From_Dictionary(input_dict, "nx_dens", ref nx_dens, true);
            Get_From_Dictionary(input_dict, "ny_dens", ref ny_dens, true);
            Get_From_Dictionary(input_dict, "nz_dens", ref nz_dens, true);
            Get_From_Dictionary(input_dict, "nx_pot", ref nx_pot, true);
            Get_From_Dictionary(input_dict, "ny_pot", ref ny_pot, true);
            Get_From_Dictionary(input_dict, "nz_pot", ref nz_pot, true);


            // and find the domain minimum coordinate values
            xmin_pot = Geom_Tool.Get_Xmin(layers);
            ymin_pot = Geom_Tool.Get_Ymin(layers);
            zmin_pot = Geom_Tool.Get_Zmin(layers);
            // but still try to get them from the dictionary if its there
            Get_From_Dictionary<double>(input_dict, "xmin_pot", ref xmin_pot, true);
            Get_From_Dictionary<double>(input_dict, "ymin_pot", ref ymin_pot, true);
            Get_From_Dictionary<double>(input_dict, "zmin_pot", ref zmin_pot, true);
            
            // and by default these are the same as for the density
            xmin_dens = xmin_pot;
            ymin_dens = ymin_pot;
            zmin_dens = zmin_pot;
            // and, once again, try to find something in the dictionary
            Get_From_Dictionary<double>(input_dict, "xmin_dens", ref xmin_dens, true);
            Get_From_Dictionary<double>(input_dict, "ymin_dens", ref ymin_dens, true);
            Get_From_Dictionary<double>(input_dict, "zmin_dens", ref zmin_dens, true);

            // and try and get as much information about the split gate dimensions as possible
            device_dimensions = new Dictionary<string, double>();
            device_dimensions.Add("top_length", Get_From_Dictionary<double>(input_dict, "top_length", 0.0));
            device_dimensions.Add("split_width", Get_From_Dictionary<double>(input_dict, "split_width", 0.0));
            device_dimensions.Add("split_length", Get_From_Dictionary<double>(input_dict, "split_length", 0.0));

            // as well as for the gate voltages
            boundary_conditions = new Dictionary<string, double>();
            boundary_conditions.Add("top_V", Get_From_Dictionary<double>(input_dict, "top_V", 0.0));
            boundary_conditions.Add("split_V1", Get_From_Dictionary<double>(input_dict, "split_V1", Get_From_Dictionary<double>(input_dict, "split_V", 0.0)));
            boundary_conditions.Add("split_V2", Get_From_Dictionary<double>(input_dict, "split_V2", Get_From_Dictionary<double>(input_dict, "split_V", 0.0)));
            boundary_conditions.Add("bottom_V", Get_From_Dictionary<double>(input_dict, "bottom_V", 0.0));
            Inputs_to_Dictionary.Parse_Voltages(input_dict, boundary_conditions);

            // work out whether we are doing dft or not
            Get_From_Dictionary<bool>(input_dict, "no_dft", ref no_dft, true);
            if (input_dict.ContainsKey("dft"))
                no_dft = !(bool)input_dict["dft"];
            // and whether to set the mixing parameter to something other than the default
            Get_From_Dictionary<double>(input_dict, "dft_mixing_parameter", ref dft_mixing_parameter, true);

            // get keys for any interesting start-up protocols
            if (input_dict.ContainsKey("hot_start")) hot_start = (bool)input_dict["hot_start"];
            Get_From_Dictionary<bool>(input_dict, "initialise_with_1D_data", ref initialise_with_1D_data, true);
            if (File.Exists("restart.flag") && input_dict.ContainsKey("with_checkpointing"))
                if (File.ReadAllLines("restart.flag")[0] == "true" && (bool)input_dict["with_checkpointing"])
                    initialise_from_restart = true;

            // and create a restart flag file if necessary
            if (!initialise_from_restart)
            {
                StreamWriter sw_flag = new StreamWriter("restart.flag"); sw_flag.WriteLine("true"); sw_flag.Close();
            }

            if (!Check_Boundary_Points(layers, Zmin_Dens, Zmin_Dens + Dz_Dens * Nz_Dens, Dz_Dens))
                throw new Exception("Error - there must be lattice points on all of the boundaries between zmin = " + Zmin_Dens.ToString() + " and zmax = " + (Zmin_Dens + Dz_Dens * Nz_Dens).ToString());

            // and load the output suffix for identification of output files
            Get_From_Dictionary<string>(input_dict, "output_suffix", ref output_suffix, true);
        }
        protected abstract void Initialise_DataClasses(Dictionary<string, object> input_dict);
        protected abstract void Initialise_from_Hot_Start(Dictionary<string, object> input_dict);

        /// <summary>
        /// The density calculations must have lattice points on all of the boundaries between zmin and zmax
        /// This method finds these boundaries and checks that the lattice has points on them
        /// </summary>
        protected bool Check_Boundary_Points(ILayer[] layers, double zmin, double zmax, double dx)
        {
            // find out how many layer iterfaces are between zmin and zmax
            int init_layer_no = Geom_Tool.GetLayer(layers, zmin).Layer_No - 1;
            int count = Geom_Tool.GetLayer(layers, zmax).Layer_No - Geom_Tool.GetLayer(layers, zmin).Layer_No;

            // check that these interfaces have lattice points on them
            for (int i = 0; i < count; i++)
                if (Math.IEEERemainder(layers[init_layer_no + i].Zmax - zmin, dx) > 1e-9)
                {
                    Console.WriteLine("WARNING - I don't think there's a lattice point on the interface at z = " + layers[init_layer_no + i].Zmax.ToString());
                    return false;
                }

            return true;
        }

        protected void Initialise_from_Restart(Dictionary<string, object> input_dict)
        {
            Console.WriteLine("Recovering data from previous, interrupted simulation...");

            // load density and chemical potential
            string[] cardens_tmp = File.ReadAllLines("carrier_density.tmp");
            string[] dopdens_tmp = File.ReadAllLines("dopent_density.tmp");
            string[] chempot_tmp = File.ReadAllLines("chem_pot.tmp");
            // and load it into the correct classes
            for (int i = 0; i < nx_dens * ny_dens * nz_dens; i++)
            {
                carrier_charge_density.Spin_Up[i] = double.Parse(cardens_tmp[i].Split(' ')[0]);
                carrier_charge_density.Spin_Down[i] = double.Parse(cardens_tmp[i].Split(' ')[1]);
                dopent_charge_density.Spin_Up[i] = double.Parse(dopdens_tmp[i].Split(' ')[0]);
                dopent_charge_density.Spin_Down[i] = double.Parse(dopdens_tmp[i].Split(' ')[1]);
                chem_pot[i] = double.Parse(chempot_tmp[i]);
            }
            // the value of t
            string[] t_tmp = File.ReadAllLines("t_val.tmp"); t = double.Parse(t_tmp[0]);

            // and load the surface charge
            boundary_conditions.Add("surface", (double)input_dict["surface_charge"]);

            Console.WriteLine("Data recovered.  Restarting from checkpoint");
        }

        public abstract bool Run();

        protected bool Run_Iteration_Routine(IDensity_Solve dens_solv, IPoisson_Solve pois_solv, double tol)
        {
            return Run_Iteration_Routine(dens_solv, pois_solv, tol, int.MaxValue);
        }
        protected abstract bool Run_Iteration_Routine(IDensity_Solve dens_solv, IPoisson_Solve pois_solv, double tol, int max_iterations);

        public void Checkpoint()
        {
            // finally, write all important data to file
            StreamWriter sw_cardens = new StreamWriter("carrier_density.tmp");
            StreamWriter sw_dopdens = new StreamWriter("dopent_density.tmp");
            StreamWriter sw_chempot = new StreamWriter("chem_pot.tmp");
            for (int i = 0; i < nx_dens * ny_dens * nz_dens; i++)
            {
                sw_cardens.WriteLine(carrier_charge_density.Spin_Up[i].ToString() + " " + carrier_charge_density.Spin_Down[i].ToString());
                sw_dopdens.WriteLine(dopent_charge_density.Spin_Up[i].ToString() + " " + dopent_charge_density.Spin_Down[i].ToString());
                sw_chempot.WriteLine(chem_pot[i].ToString());
            }
            sw_cardens.Close(); sw_dopdens.Close(); sw_chempot.Close();
            // and the current value of t
            StreamWriter sw_t = new StreamWriter("t_val.tmp"); sw_t.WriteLine(t.ToString()); sw_t.Close();
        }

        protected string Generate_Output_String(int count, Band_Data x, Band_Data dens_diff)
        {
            string output_string = "Iter = " + count.ToString() + "\tDens = " + dens_diff.Max().ToString("F4") + "\tPot = " + (Physics_Base.q_e * x.InfinityNorm()).ToString("F6") + "\tt = " + t.ToString("F6");
            return output_string;
        }

        public void Close(double unit_charge, bool converged, int no_runs)
        {
            if (!converged)
            {
                StreamWriter sw_notconverged = new StreamWriter("not_converged" + output_suffix);
                sw_notconverged.WriteLine("Not converged in " + no_runs.ToString() + " iterations");
                sw_notconverged.Close();
            }

            // save final density out
            (carrier_charge_density.Spin_Summed_Data / (unit_charge)).Save_Data("dens" + output_suffix);
            (carrier_charge_density.Spin_Up / (unit_charge)).Save_Data("dens_up" + output_suffix);
            (carrier_charge_density.Spin_Down / (unit_charge)).Save_Data("dens_down" + output_suffix);
            
            // delete the restart flag files and data
            File.Delete("restart.flag");
            File.Delete("t_val.tmp");
            File.Delete("carrier_density.tmp");
            File.Delete("dopent_density.tmp");
            File.Delete("chem_pot.tmp");
        }

        double div_fact = 0.8;
        /// <summary>
        /// calculates an optimal t based on bisection
        /// </summary>
        /// <param name="t"></param>
        /// <param name="phi"></param>
        /// <param name="x"></param>
        /// <param name="car_dens"></param>
        /// <param name="dop_dens"></param>
        /// <param name="pois_solv"></param>
        /// <param name="dens_solv"></param>
        /// <returns></returns>
        protected double Calculate_optimal_t(double t, Band_Data phi, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv, double minval)
        {
            double maxval = 3.0;
            SpinResolved_Data car_dens_copy = car_dens.DeepenThisCopy();
            SpinResolved_Data dop_dens_copy = dop_dens.DeepenThisCopy();

            double vpa = calc_vp(t, phi, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
            double vpb = calc_vp(div_fact * t, phi, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
            double t_orig = t;

            // work out whether this is going in the right direction (assuming vp is monotonic)
            if (Math.Abs(vpb) < Math.Abs(vpa))
            {
                // if 0.5 * t was going downhill, first, halve t seeing as you've already done this step
                t = div_fact * t;
                // then halve the damping parameter and check whether you've found a root yet
                while (Math.Sign(vpb) == Math.Sign(vpa))
                {
                    t = div_fact * t;
                    if (t < minval)
                        if (Math.Sign(calc_vp(1.0, phi, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv)) == Math.Sign(vpb))
                            return 0.5;
                            //return 1.0 / div_fact;
                        else
                            return Calculate_optimal_t(1.0, phi, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv, minval);
                    
                    
                    vpa = vpb;
                    vpb = calc_vp(t, phi, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                }

                //return 0.5 * (1.0 + (1.0 / div_fact)) * t;
                return t - t * (1.0 - 1.0 / div_fact) * vpb / (vpa - vpb);
            }
            else
            {
                // if 0.5 * t was going uphill, then we need to be doubling t and looking for the root
                while (Math.Sign(vpb) == Math.Sign(vpa))
                {
                    t = (1.0 / div_fact) * t;
                    if (t > maxval)
                        return maxval;

                    vpa = vpb;
                    vpb = calc_vp(t, phi, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                }

                //return 0.5 * (1.0 + div_fact) * t;
                return t - t * (1.0 - div_fact) * vpb / (vpb - vpa);
            }
             
        }

        protected virtual double calc_vp(double t, Band_Data phi, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            double vp;

            SpinResolved_Data tmp_dens = dens_solv.Get_ChargeDensity(layers, car_dens, dop_dens, Physics_Base.q_e * (phi + t * x));
            Band_Data V_Prime = -1.0 * (phi.Laplacian + t * x.Laplacian) - tmp_dens.Spin_Summed_Data;

            vp = 0.0;
            for (int i = 0; i < x.Length; i++)
            {
                vp += V_Prime[i] * x[i];
            }

            // multiply by volume element (this works in all dimensions as default dx, dy, dz are 1.0)
            return vp * dx_dens * dy_dens * dz_dens;
        }

        void Print_Data(Band_Data x, Band_Data v)
        {
            StreamWriter swx = new StreamWriter("x_prt.dat");
            StreamWriter swv = new StreamWriter("v_prt.dat");
            for (int i =0; i < v.mat.Rows; i++)
                for (int j = 0; j < v.mat.Cols; j++)
                {
                    swx.Write(x.mat[i, j].ToString() + '\t');
                    swv.Write(v.mat[i, j].ToString() + '\t');
                    if (j == v.mat.Cols - 1)
                    {
                        swv.WriteLine();
                        swx.WriteLine();
                    }
                }

            swx.Close(); swv.Close();
        }

        void Print_VP(Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            StreamWriter sw = new StreamWriter("vp");
            int count_max = 100;
            double dt = 0.01;

            for (int i = 0; i < count_max; i++)
                sw.WriteLine(calc_vp(i * dt, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv).ToString());

            sw.Close();
        }

        /// <summary>
        /// returns a list of dopent freeze-out temperatures between the initial and final temperatures of the experiment
        /// and with the final temperature at the end
        /// </summary>
        protected double[] Freeze_Out_Temperatures()
        {
            double[] raw_temps = new double[layers.Length];
            for (int i = 0; i < layers.Length; i++)
                raw_temps[i] = layers[i].Dopent_FreezeOut_T;

            // now, calculate which temperatures are in the range (final_T < T < init_T) and sort in descending order
            double[] temp_list = (from value in raw_temps where (value > temperature && value < initial_temperature) select value).ToArray().Concat(new[] {temperature}).ToArray();
            return temp_list.Distinct().OrderByDescending(c => c).ToArray();
        }

        /// <summary>
        /// gets a value of the type "T" from the input dictionary
        /// DOES NOT WORK FOR INTs
        /// </summary>
        /// <param name="attempt_only"> set to true if you don't want this method to throw an Exception</param>
        protected T Get_From_Dictionary<T>(Dictionary<string, object> input, string key, T default_value)
        {
            if (input.ContainsKey(key))
                return (T)input[key];
            else
                return default_value;
        }

        /// <summary>
        /// gets a value of the type "T" from the input dictionary
        /// DOES NOT WORK FOR INTs
        /// </summary>
        /// <param name="attempt_only"> set to true if you don't want this method to throw an Exception</param>
        protected void Get_From_Dictionary<T>(Dictionary<string, object> input, string key, ref T value, bool attempt_only = false)
        {
            if (input.ContainsKey(key))
                value = (T)input[key];
            else if (!attempt_only)
                throw new KeyNotFoundException("Error - cannot find the key: " + key);
        }

        protected void Get_From_Dictionary(Dictionary<string, object> input, string key, ref int value, bool attempt_only = false)
        {
            if (input.ContainsKey(key))
                value = (int)(double)input[key];
            else if (!attempt_only)
                throw new KeyNotFoundException("Error - cannot find the key: " + key);
        }

        protected int Get_From_Dictionary(Dictionary<string, object> input, string key, ref int value, int default_value)
        {
            if (input.ContainsKey(key))
                return (int)(double)input[key];
            else
                return default_value;
        }

        protected void Get_From_Dictionary<T>(Dictionary<string, object> input, string key, ref T value, T default_value)
        {
            Get_From_Dictionary<T>(input, key, ref default_value, true);
        }

        public SpinResolved_Data Carrier_Density
        {
            get { return carrier_charge_density; }
        }

        public SpinResolved_Data Dopent_Density
        {
            get { return dopent_charge_density; }
        }

        public Band_Data Chemical_Potential
        {
            get { return chem_pot; }
        }
        public Band_Data GPhi
        {
            get { return gphi; }
        }
        public Band_Data X
        {
            get { return x; }
        }

        public ILayer[] Layers
        {
            get { return layers; }
        }

        public double Temperature
        {
            get { return temperature; }
        }

        public double Current_Temperature
        {
            get { return current_temperature; }
        }

        public int Nx_Dens
        {
            get { return nx_dens; }
        }

        public double Dx_Dens
        {
            get { return dx_dens; }
        }

        public double Xmin_Dens
        {
            get { return xmin_dens; }
        }
        public int Ny_Dens
        {
            get { return ny_dens; }
        }

        public double Dy_Dens
        {
            get { return dy_dens; }
        }

        public double Ymin_Dens
        {
            get { return ymin_dens; }
        }
        public int Nz_Dens
        {
            get { return nz_dens; }
        }

        public double Dz_Dens
        {
            get { return dz_dens; }
        }

        public double Zmin_Dens
        {
            get { return zmin_dens; }
        }

        public int Nx_Pot
        {
            get { return nx_pot; }
        }

        public double Dx_Pot
        {
            get { return dx_pot; }
        }

        public double Xmin_Pot
        {
            get { return xmin_pot; }
        }
        public int Ny_Pot
        {
            get { return ny_pot; }
        }

        public double Dy_Pot
        {
            get { return dy_pot; }
        }

        public double Ymin_Pot
        {
            get { return ymin_pot; }
        }
        public int Nz_Pot
        {
            get { return nz_pot; }
        }

        public double Dz_Pot
        {
            get { return dz_pot; }
        }

        public double Zmin_Pot
        {
            get { return zmin_pot; }
        }
    }
}
