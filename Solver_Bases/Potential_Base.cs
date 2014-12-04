using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using CenterSpace.NMath.Core;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;
using System.Threading.Tasks;

namespace Solver_Bases
{
    public abstract class Potential_Base : IPoisson_Solve

    {
        /// this is where the density will be saved out to
        protected string dens_filename;
        protected string xc_pot_filename = "xc_pot.dat";

        protected bool flexPDE;
        protected string flexpde_location;
        protected string flexpde_inputfile;

        protected double tol;
        private bool converged;
        private double convergence_factor = double.MaxValue;

        public Potential_Base(bool using_flexPDE, string flexPDE_input, string flexPDE_location, double tol)
        {
            // check whether using flexPDE
            flexPDE = using_flexPDE;
            if (using_flexPDE)
            {
                this.flexpde_inputfile = flexPDE_input;
                this.flexpde_location = flexPDE_location;
            }

            // get the tolerance needed for the potential before saying it's good enough
            this.tol = tol;
        }

        public Band_Data Get_Chemical_Potential(Band_Data density)
        {
            if (flexPDE)
                // calculate chemical potential using a potential found by calling FlexPDE
                return Get_ChemPot_From_FlexPDE(density, dens_filename);
            else
                // calculate chemical potential using a potential calculated on a regular grid (not ideal, or scalable)
                return Get_ChemPot_On_Regular_Grid(density);
        }

        [System.Runtime.InteropServices.DllImport("user32.dll")]
        static extern int GetForegroundWindow();
        [System.Runtime.InteropServices.DllImport("user32.dll")]
        static extern int SetForegroundWindow(int handle);
        [System.Runtime.InteropServices.DllImport("user32.dll")]
        static extern int LockSetForegroundWindow(int handle);

        /// <summary>
        /// gets the band energies for the given charge distribution using flexPDE
        /// </summary>
        protected Band_Data Get_ChemPot_From_FlexPDE(Band_Data density, string dens_filename)
        {
            // save density to file in a FlexPDE "TABLE" format
            Save_Density_Data(density, dens_filename);

            // run the code
            Run_FlexPDE_Code("pot.dat");

            string[] lines = File.ReadAllLines("pot.dat");
            string[] data = Trim_Potential_File(lines);

            // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
            return Physics_Base.q_e * Parse_Potential(data);
        }

        protected virtual string[] Trim_Potential_File(string[] lines)
        {
            // work out where the data starts (this is flexPDE specific)
            int first_line = 0;
            for (int i = 0; i < lines.Length; i++)
                if (lines[i].StartsWith("}"))
                {
                    first_line = i + 1;
                    break;
                }

            // trim off the first lines which contain no data
            string[] data = new string[lines.Length - first_line];
            for (int i = first_line; i < lines.Length; i++)
                data[i - first_line] = lines[i];
            return data;
        }

        protected void Run_FlexPDE_Code(string result_filename)
        {
            // remove pot.dat if it still exists (to make sure that a new data file is made by flexPDE)
            try { File.Delete(result_filename); }
            catch (Exception) { }

            if (!File.Exists(flexpde_inputfile))
                throw new Exception("Error - there is no input file for flexpde!");

      //      Stopwatch stpwtch = new Stopwatch();
      //      stpwtch.Start();

            // run the flexPDE program as a process (quietly)
            //Process.Start("C:\\FlexPDE6\\FlexPDE6.exe", "-Q " + flexpde_inputfile);
            int handle = GetForegroundWindow();
            Process pot_process = new Process();
            pot_process.StartInfo = new ProcessStartInfo(flexpde_location, "-S " + flexpde_inputfile);
            pot_process.StartInfo.WindowStyle = ProcessWindowStyle.Minimized;

            //AutoResetEvent ev = new AutoResetEvent(false);
            //Task.Factory.StartNew(() => { pot_process.Start(); ev.Set(); });
            //while (!ev.WaitOne(0))
            //{
            //    SetForegroundWindow(handle);
            //    Thread.Sleep(0);
            //}
            int report = SetForegroundWindow(1);
            pot_process.Start();
            report = SetForegroundWindow(2);

            //Process.Start(flexpde_location, "-Q " + flexpde_inputfile);
            while (!File.Exists(result_filename))
                Thread.Sleep(100);
            Thread.Sleep(5000);

        //    stpwtch.Stop();
        //    Console.WriteLine("Time spent in FlexPDE: " + stpwtch.Elapsed.TotalSeconds.ToString() + " s");
        }

        /// <summary>
        /// Checks whether the band energies has converged by comparing old and new energies
        /// and determining whether every term is the same within a given tolerance
        /// </summary>
        public bool Check_Convergence(Band_Data band_energy, Band_Data new_band_energy, double tol)
        {
            double[] energy_diff = Get_Array_of_Absolute_Differences(band_energy, new_band_energy);

            int[] converged_test = new int[energy_diff.Length];
            for (int i = 0; i < energy_diff.Length; i++)
            {
                converged_test[i] = 0;
                if (energy_diff[i] < tol)
                    converged_test[i] = 1;
            }

            convergence_factor = energy_diff.Sum();

            if (converged_test.Sum() == energy_diff.Length)
                return true;
            else
                return false;
        }

        /// <summary>
        /// Checks whether the band energies has converged by calculating the absolute value of the blended band energies
        /// and determining whether every term is zero within a given tolerance
        /// </summary>
        public bool Check_Convergence(Band_Data blending_energy, double tol)
        {
            double[] energy_diff = new double[blending_energy.Length];
            for (int i = 0; i < blending_energy.Length; i++)
                energy_diff[i] = Math.Abs(blending_energy[i]);

            int[] converged_test = new int[energy_diff.Length];
            for (int i = 0; i < energy_diff.Length; i++)
            {
                converged_test[i] = 0;
                if (energy_diff[i] < tol)
                    converged_test[i] = 1;
            }

            convergence_factor = energy_diff.Sum();

            if (converged_test.Sum() == energy_diff.Length)
                return true;
            else
                return false;
        }

        public double Renew_Mixing_Parameter(Band_Data band_energy, Band_Data new_band_energy, double alpha_min, double alpha)
        {
            double[] energy_diff = Get_Array_of_Absolute_Differences(band_energy, new_band_energy);

            // the new mixing factor is the maximum absolute change, multiplied by the original mixing factor
            double new_mixing_parameter = alpha_min / energy_diff.Max();
            if (new_mixing_parameter > alpha_min && new_mixing_parameter < 0.01)
                return new_mixing_parameter;
            else if (new_mixing_parameter > 0.01)
                return 0.01;
            else
                return alpha_min;
        }

        /// <summary>
        /// returns an array containing |energy1 - energy2| elements
        /// </summary>
        double[] Get_Array_of_Absolute_Differences(Band_Data energy1, Band_Data energy2)
        {
            double[] result = new double[energy1.Length];
            for (int i = 0; i < energy1.Length; i++)
                result[i] = Math.Abs(energy1[i] - energy2[i]);

            return result;
        }

        public void Output(Band_Data data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            for (int i = 0; i < data.Length; i++)
                sw.WriteLine(data[i].ToString());

            sw.Close();
        }

        /// <summary>
        /// blends the old and new potentials based on a parameter using a simple mixing scheme
        /// psi_new = (1 - alpha) * psi_old + alpha * psi_calc
        /// ALSO: checks for convergence here
        /// </summary>
        public Band_Data Blend(Band_Data band_energy, Band_Data new_band_energy, double blend_parameter)
        {
            Band_Data blending_energy = band_energy - new_band_energy;

            // check for convergence
            converged = Check_Convergence(blending_energy, tol);

            return band_energy - blend_parameter * blending_energy;
        }

        /// <summary>
        /// blends the old and new band energies based on a parameter using a simple mixing scheme
        /// psi_new = (1 - alpha) * psi_old + alpha * psi_calc
        /// ALSO: checks for convergence here
        /// </summary>
        public void Blend(ref Band_Data band_energy, Band_Data new_band_energy, double blend_parameter)
        {
            Band_Data blending_energy = band_energy - new_band_energy;

            // check for convergence
            converged = Check_Convergence(blending_energy, tol);

            band_energy = band_energy - blend_parameter * blending_energy;
        }

        public bool Converged
        {
            get { return converged; }
        }

        public double Convergence_Factor
        {
            get { return convergence_factor; }
        }

        /// <summary>
        /// resets the state of the poisson solver
        /// </summary>
        public void Reset()
        {
            Console.WriteLine("Resetting convergence criteria for the potential solver");
            converged = false; convergence_factor = double.MaxValue;
        }

        public Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data rhs, SpinResolved_Data car_dens)
        {
            Save_Density_Data(car_dens.Spin_Summed_Data, dens_filename);
            Save_Density_Data(Physics_Base.Get_XC_Potential(car_dens), xc_pot_filename);
            return Calculate_Newton_Step(rho_prime, rhs);
        }

        /// <summary>
        /// temporary method!!!!
        /// </summary>
        public Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data rhs, SpinResolved_Data car_dens, Band_Data dft_diff)
        {
            Save_Density_Data(dft_diff + Physics_Base.Get_XC_Potential(car_dens), "xc_pot_calc.dat");
            return Calculate_Newton_Step(rho_prime, rhs, car_dens);
        }

        protected abstract Band_Data Parse_Potential(string[] data);
        protected abstract Band_Data Get_ChemPot_On_Regular_Grid(Band_Data density);
        protected abstract void Save_Density_Data(Band_Data density, string input_file_name);
        public abstract Band_Data Calculate_Laplacian(Band_Data input_data);
        public abstract Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data rhs);
    }
}
