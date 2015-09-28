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
        // particle-specific constants (may be changed but only through Change_Mass(double new_mass) method
        // default values are for electrons in GaAs
        protected double mass = 0.067 * Physics_Base.m_e;                 // (meV) (ps)^2 (nm)^-2 with GaAs effective mass
        protected double unit_charge = -1.0 * Physics_Base.q_e;
        protected Carrier carrier_type = Carrier.electron;

        protected double temperature;

        protected double mu = 0.0;
        protected int dim;

        private bool converged = false;
        private double convergence_factor = double.MaxValue;

        public Density_Base(double temperature)
        {
            this.temperature = temperature;

            // set the DFT parameters to the default (ie. electrons in GaAs
            Physics_Base.a_B = (4.0 * Math.PI * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar) / (mass * unit_charge * unit_charge);
            Physics_Base.Ry = (mass * unit_charge * unit_charge * unit_charge * unit_charge) / (32.0 * Math.PI * Math.PI * Physics_Base.epsilon * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar);
        }

        /// <summary>
        /// Changes the mass of the particles calculated
        /// </summary>
        public void Change_Mass(double new_mass)
        {
            Console.WriteLine("Changing the mass of the particle to new_mass = " + new_mass.ToString());
            mass = new_mass;
            Physics_Base.a_B = (4.0 * Math.PI * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar) / (mass * unit_charge * unit_charge);
            Physics_Base.Ry = (mass * unit_charge * unit_charge * unit_charge * unit_charge) / (32.0 * Math.PI * Math.PI * Physics_Base.epsilon * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar);
        }

        public void Change_Charge(double new_charge)
        {
            Console.WriteLine("Changing the charge of the particle from q_e = " + unit_charge.ToString() + " to new_charge = " + new_charge.ToString());

            // check if the new charge is +-q_e... if the relative charge is not +1 or -1, write a big message
            if (Math.Abs(new_charge) != Math.Abs(Physics_Base.q_e))
                Console.WriteLine("\n----------------------------------------------------------------\n---- The new charge is being set to a non-standard value -------\n---- This is essentially a fractionally charged system!! -------\n----------------------------------------------------------------\n");

            // change values
            unit_charge = new_charge;
            Physics_Base.a_B = (4.0 * Math.PI * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar) / (mass * unit_charge * unit_charge);
            Physics_Base.Ry = (mass * unit_charge * unit_charge * unit_charge * unit_charge) / (32.0 * Math.PI * Math.PI * Physics_Base.epsilon * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar);
        }

        /// <summary>
        /// returns the integrated density of states in the translationally invariant direction 
        /// </summary>
        protected double Get_OneD_DoS(double band_edge, double no_kb_T)
        {
            if (band_edge >= no_kb_T * Physics_Base.kB * temperature)
                return 0.0;
            else if (temperature == 0.0)
                return 2.0 * Math.Sqrt(-2.0 * mass * band_edge) / (Math.PI * Physics_Base.hbar);
            else
            {
                // calculate the density of states integral directly
                double alpha = 2.0 * Math.Sqrt(2.0 * mass) / (Math.PI * Physics_Base.hbar * Physics_Base.kB * temperature);
                double beta = 1.0 / (Physics_Base.kB * temperature);
                OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double E) => Math.Sqrt(E - band_edge) * Math.Exp(beta * E) * Math.Pow(Math.Exp(beta * E) + 1, -2.0)));
                return Perform_DoS_Integral(band_edge, no_kb_T, alpha, dos_integrand);
            }
        }

        protected double Get_OneD_DoS_Deriv(double band_edge, double no_kb_T)
        {
            if (band_edge >= no_kb_T * Physics_Base.kB * temperature)
                return 0.0;
            else if (temperature == 0.0)
                return -1.0 * Math.Sqrt(-2.0 * mass / band_edge) / (Math.PI * Physics_Base.hbar);
            else
            {
                // calculate the density of states integral directly
                double alpha = 2.0 * Math.Sqrt(2.0 * mass) / (Math.PI * Physics_Base.hbar * Physics_Base.kB * temperature);
                double beta = 1.0 / (Physics_Base.kB * temperature);
                OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double E) => Math.Sqrt(E - band_edge) * Math.Exp(beta * E) * Math.Pow(Math.Exp(beta * E) + 1, -3.0) * (beta * (Math.Exp(beta * E) + 1) - 2.0 * beta * Math.Exp(beta * E))));
                return Perform_DoS_Integral(band_edge, no_kb_T, alpha, dos_integrand);
            }
        }

        protected double Get_TwoD_DoS(double band_edge, double no_kb_T)
        {
            if (band_edge >= no_kb_T * Physics_Base.kB * temperature)
                return 0.0;
            else if (temperature == 0.0)
                return -1.0 * band_edge * mass / (Physics_Base.hbar * Physics_Base.hbar * 2.0 * Math.PI);
            else
            {
                // calculate the density of states integral directly
                double alpha = mass / (Physics_Base.hbar * Physics_Base.hbar * 2.0 * Math.PI);
                double beta = 1.0 / (Physics_Base.kB * temperature);
                OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double E) => 1.0 / (Math.Exp(beta * E) + 1)));
                return Perform_DoS_Integral(band_edge, no_kb_T, alpha, dos_integrand);
            }
        }

        protected double Get_TwoD_DoS_Deriv(double band_edge, double no_kb_T)
        {
            if (band_edge >= no_kb_T * Physics_Base.kB * temperature)
                return 0.0;
            else if (temperature == 0.0)
                return mass / (Physics_Base.hbar * Physics_Base.hbar * 2.0 * Math.PI);
            else
            {
                // calculate the derivative of the density of states integral directly
                // NOTE: The 2D DoS is constant so this is just the integral of the derivative of the Fermi function
                double alpha = mass / (Physics_Base.hbar * Physics_Base.hbar * 2.0 * Math.PI);
                double beta = 1.0 / (Physics_Base.kB * temperature);
                OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double E) => beta * Math.Exp(beta * E) * Math.Pow(Math.Exp(beta * E) + 1, -2.0)));
                return Perform_DoS_Integral(band_edge, no_kb_T, alpha, dos_integrand);
            }
        }

        private double Perform_DoS_Integral(double band_edge, double no_kb_T, double alpha, OneVariableFunction dos_integrand)
        {
            dos_integrand.Integrator = new GaussKronrodIntegrator(); 

            if (band_edge < -1.0 * no_kb_T * Physics_Base.kB * temperature)
                return alpha * dos_integrand.Integrate(-1.0 * no_kb_T * Physics_Base.kB * temperature, no_kb_T * Physics_Base.kB * temperature);
            else
                return alpha * dos_integrand.Integrate(band_edge, no_kb_T * Physics_Base.kB * temperature);
        }

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

        #region fermifunctions
        /// <summary>
        /// Gets the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function(double energy)
        {
            return Get_Fermi_Function(energy, mu, temperature);
        }

        /// <summary>
        /// Calculates the fermi function for arbitrary energy, E_f and T
        /// </summary>
        protected double Get_Fermi_Function(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy > mu)
                    return 0.0;
                else
                    return 1.0;
            else
                return 1.0 / (Math.Exp((energy - mu) / (Physics_Base.kB * T)) + 1.0);
        }

        /// <summary>
        /// Gets the derivative of the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function_Derivative(double energy)
        {
            return Get_Fermi_Function_Derivative(energy, mu, temperature);
        }

        /// <summary>
        /// Calculates the derivative of the fermi function for arbitrary energy, E_f and T
        /// </summary>
        protected double Get_Fermi_Function_Derivative(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy == mu)
                    return 1.0;
                else
                    return 0.0;
            else
            {
                double exponent = (energy - mu) / (Physics_Base.kB * T);

                if (double.IsInfinity(Math.Exp(exponent)))
                    return 0.0;
                else
                    return Math.Exp(exponent) / (Math.Pow((Math.Exp(exponent) + 1.0), 2.0) * Physics_Base.kB * T);
            }
        }

        /// <summary>
        /// Gets the fermi function for the dopents at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Dopent_Fermi_Function(double energy)
        {
            return Get_Dopent_Fermi_Function(energy, mu, temperature);
        }

        /// <summary>
        /// Calculates the spin-resolved fermi function for a dopent.
        /// This is different from the typical fermi function as double occupation of the donor is not allowed
        /// </summary>
        protected double Get_Dopent_Fermi_Function(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy > mu)
                    return 0.0;
                else
                    return 1.0;
            else
                return 1.0 / (Math.Exp((energy - mu) / (Physics_Base.kB * T)) + 2.0);
        }

        /// <summary>
        /// Gets the derivative of the fermi function for the dopents at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Dopent_Fermi_Function_Derivative(double energy)
        {
            return Get_Dopent_Fermi_Function_Derivative(energy, mu, temperature);
        }

        /// <summary>
        /// Calculates the derivative of the spin-resolved fermi function for a dopent.
        /// This is different from the typical fermi function as double occupation of the donor is not allowed
        /// </summary>
        protected double Get_Dopent_Fermi_Function_Derivative(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy == mu)
                    return 1.0;
                else
                    return 0.0;
            else
            {
                double exponent = (energy - mu) / (Physics_Base.kB * T);

                if (double.IsInfinity(Math.Exp(exponent)))
                    return 0.0;
                else
                    return Math.Exp(exponent) / (Math.Pow((Math.Exp(exponent) + 2.0), 2.0) * Physics_Base.kB * T);
            }
        }
        #endregion

        /// <summary>
        /// for this method, it is assumed that the dopent density is frozen out (and therefore not altered) unless overriden
        /// </summary>
        public virtual void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_charge_density, ref SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            Get_ChargeDensity(layers, ref carrier_charge_density, chem_pot);
        }

        public bool Converged
        {
            get { return converged; }
        }

        public double Convergence_Factor
        {
            get { return convergence_factor; }
        }

        public abstract void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot);
        public abstract SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot);
        public abstract SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot);
        public abstract DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot);

        public void Close()
        {
            Console.WriteLine("Closing density solver");
        }

        #region dftpot
        protected Band_Data dft_pot;
        protected double alpha_dft = 0.1;
        public void Set_DFT_Potential(SpinResolved_Data car_dens)
        {
            if (this.dft_pot == null && alpha_dft != 0.0)
                this.dft_pot = Physics_Base.Get_XC_Potential(car_dens);
            else if (alpha_dft == 0.0)
                this.dft_pot = 0.0 * car_dens.Spin_Summed_Data.DeepenThisCopy();
            else
                this.dft_pot = (1.0 - alpha_dft) * dft_pot + alpha_dft * Physics_Base.Get_XC_Potential(car_dens);
        }
        public void Reset_DFT_Potential()
        {
            this.dft_pot = null;
        }
        public void Print_DFT_diff(SpinResolved_Data car_dens)
        {
            Console.WriteLine("Maximum absolute difference in DFT potentials = " + DFT_diff(car_dens).InfinityNorm().ToString("F6"));
        }
        public Band_Data DFT_diff(SpinResolved_Data car_dens)
        {
            if (alpha_dft == 0.0)
                return 0.0 * dft_pot;
            else
                return dft_pot - Physics_Base.Get_XC_Potential(car_dens);
        }
        public double DFT_Mixing_Parameter
        {
            get { return alpha_dft; }
            set { alpha_dft = value; }
        }
        #endregion
    }
}
