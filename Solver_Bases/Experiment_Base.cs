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
    public abstract class Experiment_Base
    {
        protected SpinResolved_Data carrier_density;
        protected SpinResolved_Data dopent_density;
        protected SpinResolved_Data carrier_density_deriv;
        protected SpinResolved_Data dopent_density_deriv;
        protected Band_Data chem_pot;
        protected ILayer[] layers;

        // damping parameter for Bank-Rose convergence method
        protected double t = 0.5;

        protected double bottom_V;

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
        protected string flexPDE_input;
        protected string flexPDE_location;

        protected double initial_temperature = 300.0;
        protected double temperature;

        protected bool TF_only = false;       // only run Thomas-Fermi density calculations (i.e. no DFT, Green's functions etc...)

        protected void Initialise(Dictionary<string, object> input_dict)
        {
            // solver inputs
            Get_From_Dictionary<double>(input_dict, "tolerance", ref tol);
            Get_From_Dictionary<double>(input_dict, "alpha", ref alpha);
            Get_From_Dictionary<double>(input_dict, "alpha", ref alpha_prime);

            // will not use FlexPDE unless told to
            if (input_dict.ContainsKey("use_FlexPDE")) this.using_flexPDE = (bool)input_dict["use_FlexPDE"]; else using_flexPDE = false;
            // default input file for FlexPDE is called "default.pde"
            if (input_dict.ContainsKey("FlexPDE_file")) this.flexPDE_input = (string)input_dict["FlexPDE_file"]; else this.flexPDE_input = "default.pde";
            if (using_flexPDE)
            {
                // make sure that FlexPDE does exist at the specified location
                try { this.flexPDE_location = (string)input_dict["FlexPDE_location"]; }
                catch (Exception) { throw new Exception("Error - no location for FlexPDE executable supplied"); }
                if (!File.Exists(flexPDE_location))
                    throw new Exception("Error - FlexPDE executable file does not exist at location" + flexPDE_location + "!");
            }

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
        }

        public abstract void Run();
        /*
        protected double Calculate_optimal_t(double t, Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            double vp = calc_vp(t, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);

            // search for best value of t
            if (vp > 0.0)
                while (calc_vp(2.0 * t, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv) > 0.0)
                    t = 2.0 * t;
            else
                while (calc_vp(t, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv) < 0.0)
                    t = 0.5 * t;

            return t;
        }*/

        /* 
        /// <summary>
        /// calculates the optimal t using a Brent-Dekker method (see Wikipedia)
        /// </summary>
        /// 
        protected double Calculate_optimal_t(double t, Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            double a = 0.0; double b = 1.0;
            double delta = 0.001;
            int max_iteration = 100; int count = 0;

            double vpa = calc_vp(a, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);
            double vpb = calc_vp(b, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);

            if (Math.Abs(vpa) < Math.Abs(vpb))
            {
                double tmp = a;
                a = b;
                b = tmp;
                tmp = vpa;
                vpa = vpb;
                vpb = tmp;
            }

            double c = a;
            double vpc = vpa;

            double s = 0.0; bool flag = true; double d = 2.0;
            while (!(vpb == 0) && Math.Abs(a - b) > delta && max_iteration > count)
            {
                if (vpa != vpc && vpb != vpc)
                    s = (a * vpb * vpc) / (vpa - vpb) / (vpa - vpc) + (b * vpa * vpc) / (vpb - vpa) / (vpb - vpc) + (c * vpa * vpb) / (vpc - vpa) / (vpc - vpb);
                else
                    s = b - vpb * (b - a) / (vpb - vpa);

                if (!((s > (3 * a + b) / 4 && s < b) || ((s < (3 * a + b) / 4 && s > b)))
                    || (flag && Math.Abs(s - b) >= Math.Abs(b - c) / 2)
                    || (!flag && Math.Abs(s - b) >= Math.Abs(c - d) / 2))
                {
                    s = (a + b) / 2;
                    flag = true;
                }
                else
                {
                    if ((flag && Math.Abs(b - c) < delta) || (!flag && Math.Abs(c - d) < delta))
                    {
                        s = (a + b) / 2;
                        flag = true;
                    }
                    else
                        flag = false;
                }

                double vps = calc_vp(s, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);

                d = c;
                c = b;
                vpc = vpb;

                if (vpa * vps < 0)
                {
                    b = s; vpb = vps;
                }
                else
                {
                    a = s; vpa = vps;
                }

                if (Math.Abs(vpa) < Math.Abs(vpb))
                {
                    double tmp = a;
                    a = b;
                    b = tmp;
                    tmp = vpa;
                    vpa = vpb;
                    vpb = tmp;
                }

                count++;
            }

            return b;
        }
         */

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
            double vpb = calc_vp(0.5 * t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);

            // work out whether this is going in the right direction (assuming vp is monotonic)
            if (Math.Abs(vpb) < Math.Abs(vpa))
            {
                // if 0.5 * t was going downhill, first, halve t seeing as you've already done this step
                t = 0.5 * t;
                // then halve the damping parameter and check whether you've found a root yet
                while (Math.Sign(vpb) == Math.Sign(vpa))
                {
                    if (t < minval)
                        return minval;
                    t = 0.5 * t;
                    vpa = vpb;
                    vpb = calc_vp(t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                }

                return 1.5 * t;
            }
            else
            {
                // if 0.5 * t was going downhill, then we need to be doubling t and looking for the root
                while (Math.Sign(vpb) == Math.Sign(vpa) && t < maxval)
                {
                    if (t > maxval)
                        return maxval;
                    t = 2.0 * t;
                    vpa = vpb;
                    vpb = calc_vp(2.0 * t, band_energy, x, car_dens_copy, dop_dens_copy, pois_solv, dens_solv);
                }

                return 0.75 * t;
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

        public double Bottom_BC
        {
            get { return bottom_V; }
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
