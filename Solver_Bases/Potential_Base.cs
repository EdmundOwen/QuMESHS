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

namespace Solver_Bases
{
    public abstract class Potential_Base
    {
        /// this is where the density will be saved out to
        protected string dens_filename;

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

        /// <summary>
        /// gets the band energies for the given charge distribution using flexPDE
        /// </summary>
        protected Band_Data Get_ChemPot_From_FlexPDE(Band_Data density, string dens_filename)
        {
            // save density to file in a FlexPDE "TABLE" format
            Save_Density_Data(density, dens_filename);

            // remove pot.dat if it still exists (to make sure that a new data file is made by flexPDE)
            try { File.Delete("pot.dat"); }
            catch (Exception) { }

            if (!File.Exists(flexpde_inputfile))
                throw new Exception("Error - there is no input file for flexpde!");

            // run the flexPDE program as a process (quietly)
            //Process.Start("C:\\FlexPDE6\\FlexPDE6.exe", "-Q " + flexpde_inputfile);
            Process.Start(flexpde_location, "-Q " + flexpde_inputfile);
            while (!File.Exists("pot.dat"))
                Thread.Sleep(10);
            Thread.Sleep(1000);

            string[] lines = File.ReadAllLines("pot.dat");

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

            // and trim all of the empty lines
            //string[] data = (from items in tmp where items != "" select items).ToArray();

            // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
            return Physics_Base.q_e * Parse_Potential(data);
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
        /// blends the old and new potentials based on a parameter using a simple (I think it's a Newton) scheme
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
        /// blends the old and new band energies based on a parameter using a simple (I think it's a Newton) scheme
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

        protected abstract Band_Data Parse_Potential(string[] data);
        protected abstract Band_Data Get_ChemPot_On_Regular_Grid(Band_Data density);
        protected abstract void Save_Density_Data(Band_Data density, string input_file_name);
    }
}
