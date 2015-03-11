﻿using System;
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
        protected SpinResolved_Data carrier_density;
        protected SpinResolved_Data dopent_density;
        protected SpinResolved_Data carrier_density_deriv;
        protected SpinResolved_Data dopent_density_deriv;
        protected Band_Data chem_pot;
        protected ILayer[] layers;

        // damping parameter for Bank-Rose convergence method
        protected double t = 0.5;

        // parameters for the density domain
        protected double dx_dens = 1.0, dy_dens = 1.0, dz_dens = 1.0;
        protected double xmin_dens, ymin_dens, zmin_dens = -1.0;
        protected int nx_dens, ny_dens, nz_dens = 1;

        // parameters for the potential domain
        protected double dx_pot, dy_pot, dz_pot;
        protected double xmin_pot, ymin_pot, zmin_pot = -1.0;
        protected int nx_pot, ny_pot, nz_pot;

        protected double alpha, alpha_prime, tol;

        protected bool using_flexPDE = false;
        protected bool using_dealii = false;

        protected double initial_temperature = 300.0;
        protected double temperature;

        protected bool no_dft = false;       // do not run with dft potential (ie. Hartree approximation)
        protected bool hot_start = false;     // am the program starting from a precalculated density or do i start from scratch

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // solver inputs
            Get_From_Dictionary<double>(input_dict, "tolerance", ref tol);
            Get_From_Dictionary<double>(input_dict, "alpha", ref alpha);
            Get_From_Dictionary<double>(input_dict, "alpha", ref alpha_prime);

            // will not use FlexPDE unless told to
            if (input_dict.ContainsKey("use_FlexPDE")) this.using_flexPDE = (bool)input_dict["use_FlexPDE"]; else using_flexPDE = false;
            // and the same for dealii
            if (input_dict.ContainsKey("use_deal.II")) this.using_dealii = (bool)input_dict["use_deal.II"]; else using_dealii = false;

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

            // work out whether we are doing dft or not
            // first check that there is only one entry (else throw exception)
            if ((input_dict.ContainsKey("no_dft") && input_dict.ContainsKey("dft")) || (input_dict.ContainsKey("no_dft") && input_dict.ContainsKey("TF_only")) || (input_dict.ContainsKey("dft") && input_dict.ContainsKey("TF_only")))
                throw new Exception("Error - more that one key specifying whether dft should be used!");

            // then input from dictionary
            Get_From_Dictionary<bool>(input_dict, "TF_only", ref no_dft, true);
            Get_From_Dictionary<bool>(input_dict, "no_dft", ref no_dft, true);
            if (input_dict.ContainsKey("dft"))
                no_dft = !(bool)input_dict["dft"];
        }

        public abstract void Run();

        double div_fact = 0.8;
        /// <summary>
        /// calculates an optimal t based on bisection
        /// </summary>
        /// <param name="t"></param>
        /// <param name="band_energy"></param>
        /// <param name="x"></param>
        /// <param name="car_dens"></param>
        /// <param name="dop_dens"></param>
        /// <param name="pois_solv"></param>
        /// <param name="dens_solv"></param>
        /// <returns></returns>
        protected double Calculate_optimal_t(double t, Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv, double minval)
        {
            double maxval = 1.0;
            SpinResolved_Data car_dens_copy = car_dens.DeepenThisCopy();
            SpinResolved_Data dop_dens_copy = dop_dens.DeepenThisCopy();

            double vpa = calc_vp(t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
            double vpb = calc_vp(div_fact * t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
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
                    {
                        // evaluate vp at zero
                        double vpc = calc_vp(0.0, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                        // if there is a root between minval and zero, return minval
                        if (Math.Sign(vpb) != Math.Sign(vpc))
                            return minval;
                        /*
                    // otherwise, return a negative value which can be used as a flag or, if not, will still
                    // generate weird behavious which should push the iterator into an active regime
                    else
                        return -1.0 * minval;
                         */
                        // new behaviour is to just set the blend parameter to the minimum value so that at least something happens
                        else
                            return minval;
                    }

                    vpa = vpb;
                    vpb = calc_vp(t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                }

                // check that t is not less than the minimum possible value (just in case)
                if (t < minval)
                    return minval;

                return 0.5 * ((1.0 + div_fact) / div_fact) * t;
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
                    vpb = calc_vp((1.0 / div_fact) * t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                }

                return 0.5 * (1.0 + div_fact) * t;
            }
        }

        protected virtual double calc_vp(double t, Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            double vp;

            SpinResolved_Data tmp_dens = dens_solv.Get_ChargeDensity(layers, car_dens, dop_dens, band_energy + t * x);
            Band_Data V_Prime = -1.0 * pois_solv.Calculate_Laplacian((band_energy + t * x) / Physics_Base.q_e) - tmp_dens.Spin_Summed_Data;

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
            double dt = 0.001;

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

        protected void Get_From_Dictionary<T>(Dictionary<string, object> input, string key, ref T value, T default_value)
        {
            Get_From_Dictionary<T>(input, key, ref default_value, true);
        }

        public SpinResolved_Data Carrier_Density
        {
            get { return carrier_density; }
        }

        public SpinResolved_Data Dopent_Density
        {
            get { return dopent_density; }
        }

        public Band_Data Chemical_Potential
        {
            get { return chem_pot; }
        }

        public ILayer[] Layers
        {
            get { return layers; }
        }

        public double Temperature
        {
            get { return temperature; }
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
