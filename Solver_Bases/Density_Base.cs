using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public abstract class Density_Base : IDensity_Solve
    {
        protected double temperature;

        protected double dx, dy, dz;
        protected double xmin, ymin, zmin;
        protected int nx, ny, nz;

        protected double fermi_Energy = 0.0;
        protected int dim;

        private bool converged = false;
        private double convergence_factor = double.MaxValue;

        public Density_Base(double temperature, double dx, double dy, double dz, int nx, int ny, int nz, double xmin, double ymin, double zmin)
        {
            this.temperature = temperature;
            this.dx = dx; this.dy = dy; this.dz = dz;
            this.nx = nx; this.ny = ny; this.nz = nz;
            this.xmin = xmin; this.ymin = ymin; this.zmin = zmin;

            // check how many actual dimensions are present
            if (nx != 1)
                dim = 3;
            else if (ny != 1)
                dim = 2;
            else
                dim = 1;
        }

        /// <summary>
        /// Gets the dopent occupation probability at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_NonSpinResolved_Dopent_Occupation(double energy)
        {
            return Get_NonSpinResolved_Dopent_Occupation(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Calculates the dopent occupation probability for arbitrary energy, E_f and T
        /// </summary>
        protected double Get_NonSpinResolved_Dopent_Occupation(double energy, double E_f, double T)
        {
            if (T == 0)
                if (energy > E_f)
                    return 0.0;
                else
                    return 1.0;
            else
                return 2.0 / (Math.Exp((energy - fermi_Energy) / (Physics_Base.kB * T)) + 2.0);
        }

        /// <summary>
        /// returns the integrated density of states in the translationally invariant direction 
        /// </summary>
        protected double Get_OneD_DoS(double band_edge, double no_kb_T)
        {
            if (band_edge > no_kb_T * Physics_Base.kB * temperature)
                return 0.0;
            else if (temperature == 0)
                return 2.0 * Math.Sqrt(-2.0 * Physics_Base.mass * band_edge) / (Math.PI * Physics_Base.hbar);
            else
            {
                // calculate the density of states integral directly
                double alpha = 2.0 * Math.Sqrt(2.0 * Physics_Base.mass) / (Math.PI * Physics_Base.hbar * Physics_Base.kB * temperature);
                double beta = 1.0 / (Physics_Base.kB * temperature);
                OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double E) => Math.Sqrt(E - band_edge) * Math.Exp(beta * E) * Math.Pow(Math.Exp(beta * E) + 1, -2.0)));
                dos_integrand.Integrator = new GaussKronrodIntegrator();
                return alpha * dos_integrand.Integrate(band_edge, no_kb_T * Physics_Base.kB * temperature);
            }
        }

        /*public void Output(SpinResolved_DoubleVector data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            // output the charge density
            for (int i = 0; i < data.Nx; i++)
                sw.WriteLine(data.Spin_Summed_Vector[i].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_DoubleMatrix data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - Ordering compared to Band_Data objects is not guaranteed!");
            
            // output the charge density
            for (int i = 0; i < data.Nx; i++)
                for (int j = 0; j < data.Ny; j++)
                    sw.WriteLine(data.Spin_Summed_Matrix[i, j].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_DoubleMatrix[] data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - Ordering compared to Band_Data objects is not guaranteed!");

            // output the charge density
            for (int i = 0; i < data[0].Nx; i++)
                for (int j = 0; j < data[0].Ny; j++)
                    for (int k = 0; k < data.Length; k++)
                        sw.WriteLine(data[k].Spin_Summed_Matrix[i, j].ToString());

            sw.Close();
        }*/

        public void Output(SpinResolved_Data data, string filename, bool with_warnings)
        {
            StreamWriter sw = new StreamWriter(filename);
            if (with_warnings)
            {
                sw.WriteLine("Warning - The data has been written out serially and there is no information as to which order the dimensions come in.");
                sw.WriteLine("Warning - Ordering compared to Band_Data objects is not guaranteed!");
            }

            // output the charge density
            Band_Data tot_charge = data.Spin_Summed_Data;
            for (int i = 0; i < tot_charge.Length; i++)
                sw.WriteLine(tot_charge[i].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_Data data, string filename)
        { Output(data, filename, true); }

        /// <summary>
        /// Gets the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function(double energy)
        {
            return Physics_Base.Get_Fermi_Function(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Gets the fermi function for the dopents at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Dopent_Fermi_Function(double energy)
        {
            return Physics_Base.Get_Dopent_Fermi_Function(energy, fermi_Energy, temperature);
        }

        public virtual double Get_Chemical_Potential(double z, ILayer[] layers)
        {
            return Get_Chemical_Potential(z, layers, temperature);
        }

        public virtual double Get_Chemical_Potential(double y, double z, ILayer[] layers)
        {
            return Get_Chemical_Potential(y, z, layers, temperature);
        }

        public virtual double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers)
        {
            return Get_Chemical_Potential(x, y, z, layers, temperature);
        }

        /// <summary>
        /// set the default chemical potential to be calculated at (0,0,z).  Everything else requires this method to be overridden
        /// </summary>
        public virtual double Get_Chemical_Potential(double z, ILayer[] layers, double temperature_input)
        {
            return Get_Chemical_Potential(0.0, 0.0, z, layers, temperature_input);
        }

        /// <summary>
        /// set the default chemical potential to be calculated at (0,y,z).  Everything else requires this method to be overridden
        /// </summary>
        public virtual double Get_Chemical_Potential(double y, double z, ILayer[] layers, double temperature_input)
        {
            return Get_Chemical_Potential(0.0, y, z, layers, temperature_input);
        }

        /// <summary>
        /// returns an array containing |density1 - density2| elements
        /// </summary>
        double[] Get_Array_of_Absolute_Differences(SpinResolved_Data density1, SpinResolved_Data density2)
        {
            double[] result = new double[density1.Length];
            for (int i = 0; i < density1.Spin_Summed_Data.Length; i++)
                result[i] = Math.Abs(density1.Spin_Summed_Data[i] - density2.Spin_Summed_Data[i]);

            return result;
        }

        /// <summary>
        /// Checks whether the density has converged by calculating the absolute value of the blended densitiies
        /// and determining whether every term is zero within a given tolerance
        /// </summary>
        public bool Check_Convergence(SpinResolved_Data blending_density, double tol)
        {
            double[] density_diff = new double[blending_density.Spin_Summed_Data.Length];
            for (int i = 0; i < blending_density.Spin_Summed_Data.Length; i++)
                density_diff[i] = Math.Abs(blending_density.Spin_Summed_Data[i]);

            int[] converged_test = new int[density_diff.Length];
            for (int i = 0; i < density_diff.Length; i++)
            {
                converged_test[i] = 0;
                if (density_diff[i] < tol)
                    converged_test[i] = 1;
            }

            convergence_factor = density_diff.Sum();

            if (converged_test.Sum() == density_diff.Length)
                return true;
            else
                return false;
        }

        /// <summary>
        /// Checks whether the density has converged by calculating the fractional absolute value of the blended densities
        /// and determining whether every term is zero within a given tolerance
        /// </summary>
        public bool Check_Convergence_Fraction(SpinResolved_Data blending_density, SpinResolved_Data density, double tol)
        {
            // minimum density where convergence is no longer checked is equal to maximum fluctuation
            double dens_max = Math.Abs(density.Spin_Summed_Data.mat.Min());
            double dens_min = Math.Abs(density.Spin_Summed_Data.mat.Max());
            double minval = Math.Abs(Math.Max(dens_min, dens_max) * tol);

            double[] density_diff = new double[blending_density.Spin_Summed_Data.Length];
            for (int i = 0; i < blending_density.Spin_Summed_Data.Length; i++)
                if (Math.Abs(density.Spin_Summed_Data[i]) > minval)
                    density_diff[i] = Math.Abs(blending_density.Spin_Summed_Data[i]) / Math.Abs(density.Spin_Summed_Data[i]);
                else
                    density_diff[i] = 0;

            int[] converged_test = new int[density_diff.Length];
            for (int i = 0; i < density_diff.Length; i++)
            {
                converged_test[i] = 0;
                if (density_diff[i] < tol)
                    converged_test[i] = 1;
            }

            // the convergence factor is the sum of the absolute values of the density error defined by the
            // difference between the current density and the ideal one for this potential
            convergence_factor = blending_density.Spin_Summed_Data.mat.Select(x => Math.Abs(x)).ToList().Sum();

            if (converged_test.Sum() == density_diff.Length)
                return true;
            else
                return false;
        }

        /// <summary>
        /// blends the old and new densities based on a parameter using a simple (I think it's a Newton) scheme
        /// rho_new = (1 - alpha) * rho_old + alpha * rho_calc
        /// ALSO: checks for convergence here
        /// </summary>
        public void Blend(ref SpinResolved_Data band_density, SpinResolved_Data new_band_density, double blend_parameter, double tol)
        {
            SpinResolved_Data blending_density = band_density - new_band_density;

            // check for convergence
            converged = Check_Convergence(blending_density, tol);

            band_density = band_density - blend_parameter * blending_density;
        }

        /// <summary>
        /// uses a three step blending to calculate a "damped" mixture of the densities... see notes
        /// </summary>
        public void Blend(ref SpinResolved_Data band_density, ref SpinResolved_Data old_band_density, SpinResolved_Data new_band_density, double alpha, double zeta, double tol)
        {
            SpinResolved_Data tmp_band_density = band_density.DeepenThisCopy();
            //band_density = ((2 - alpha * alpha) * tmp_band_density - (1 - zeta * alpha) * old_band_density + (alpha * alpha) * new_band_density) / (1 + zeta * alpha);
            band_density = (2 - zeta) * tmp_band_density - (1 - zeta) * old_band_density + (alpha * alpha) * (new_band_density - tmp_band_density);

            // check for convergence
            //converged = Check_Convergence_Fraction(new_band_density - band_density, band_density, tol);

            // absolute tolerance is the fractional tolerance (inputted) at the maximum density
            double dens_max = Math.Abs(band_density.Spin_Summed_Data.mat.Min());
            double dens_min = Math.Abs(band_density.Spin_Summed_Data.mat.Max());
            double tmp_tol = Math.Abs(Math.Max(dens_min, dens_max) * tol);
            converged = Check_Convergence(new_band_density - band_density, tmp_tol);

            // and reset the "old_band_density" to the new previous time step
            old_band_density = tmp_band_density;
        }

        public bool Converged
        {
            get { return converged; }
        }

        public double Convergence_Factor
        {
            get { return convergence_factor; }
        }

        public abstract void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
        public abstract void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_density, ref SpinResolved_Data dopent_density, Band_Data chem_pot);
        public abstract SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot);
        public abstract double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input);
        public abstract void Close();
    }
}
