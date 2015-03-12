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
        protected string xc_pot_filename = "xc_pot.dat";

        protected bool external_code;
        protected string external_arguments = "";

        protected string initcalc_location;
        protected string newton_location;

        protected string initcalc_result_filename = "pot.dat";
        protected string newton_result_filename = "x.dat";

        private bool converged;
        private double convergence_factor = double.MaxValue;

        double t;

        public Potential_Base(bool using_external_code)
        {
            // check whether using an external code
            external_code = using_external_code;
        }

        public Band_Data Get_Chemical_Potential(Band_Data density)
        {
            if (external_code)
                // calculate chemical potential using a potential found by calling external code
                return Get_ChemPot_From_External(density);
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
        /// gets the band energies for the given charge distribution using external code
        /// </summary>
        protected Band_Data Get_Data_From_External(string external_program_location, string external_program_argument, string result_filename)
        {
            // run the code
            Run_External_Code(external_program_location, external_program_argument, result_filename);

            string[] lines = File.ReadAllLines(result_filename);
            string[] data = Trim_Potential_File(lines);

            // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
            return Physics_Base.q_e * Parse_Potential(data);
        }

        /// <summary>
        /// run external code with null arguments, waiting until the result_filename file is created
        /// </summary>
        protected Band_Data Get_Data_From_External(string external_program_location, string result_filename)
        {
            return Get_Data_From_External(external_program_location, "", result_filename);
        }

        /// <summary>
        /// run external code with input arguments, waiting until the result_filename file is created
        /// </summary>
        protected void Run_External_Code(string external_program_location, string external_program_arguments, string result_filename)
        {
            // remove result file if it still exists (to make sure that a new data file is made by the external program)
            try { File.Delete(result_filename); }
            catch (Exception) { }

            // check whether the program is actually there
            if (!File.Exists(external_program_location))
                throw new Exception("Error - there is no file at the requested location: " + external_program_location);

            Process pot_process = new Process();
            // input arguments into process
            pot_process.StartInfo = new ProcessStartInfo(external_program_location, external_program_arguments);

            // and run
            pot_process.Start();

            // wait until the requested file is available, then exit
            while (!File.Exists(result_filename))
                Thread.Sleep(100);
            Thread.Sleep(5000);
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

        public double T { get { return t; } set { t = value; } }

        public virtual Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data carrier_density, Band_Data dft_difference)
        {
            throw new NotImplementedException();
        }

        public abstract void Initiate_Poisson_Solver(Dictionary<string, double> device_dimensions, Dictionary<string, double> boundary_conditions);
        protected abstract string[] Trim_Potential_File(string[] lines);
        protected abstract Band_Data Parse_Potential(string[] data);
        protected abstract Band_Data Get_ChemPot_From_External(Band_Data density);
        protected abstract Band_Data Get_ChemPot_On_Regular_Grid(Band_Data density);
        public abstract Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi);
        protected abstract void Save_Data(Band_Data density, string input_file_name);
        public abstract Band_Data Calculate_Laplacian(Band_Data input_data);
        public abstract Band_Data Chemical_Potential { get; }
    }
}
