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
        protected double a_B;       // Bohr radius for GaAs in nm
        protected double Ry;     // rydberg in meV in GaAs

        protected double temperature;

        protected double mu = 0.0;
        protected int dim;

        private bool converged = false;
        private double convergence_factor = double.MaxValue;

        public Density_Base(double temperature)
        {
            this.temperature = temperature;

            // set the DFT parameters to the default (ie. electrons in GaAs
            a_B = (4.0 * Math.PI * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar) / (mass * unit_charge * unit_charge);
            Ry = (mass * unit_charge * unit_charge * unit_charge * unit_charge) / (32.0 * Math.PI * Math.PI * Physics_Base.epsilon * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar);
        }

        /// <summary>
        /// Changes the mass of the particles calculated
        /// </summary>
        public void Change_Mass(double new_mass)
        {
            Console.WriteLine("Changing the mass of the particle to new_mass = " + new_mass.ToString());
            mass = new_mass;
            a_B = (4.0 * Math.PI * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar) / (mass * unit_charge * unit_charge);
            Ry = (mass * unit_charge * unit_charge * unit_charge * unit_charge) / (32.0 * Math.PI * Math.PI * Physics_Base.epsilon * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar);
        }

        public void Change_Charge(double new_charge)
        {
            Console.WriteLine("Changing the charge of the particle from q_e = " + unit_charge.ToString() + " to new_charge = " + new_charge.ToString());

            // check if the new charge is +-q_e... if the relative charge is not +1 or -1, write a big message
            if (Math.Abs(new_charge) != Math.Abs(Physics_Base.q_e))
                Console.WriteLine("\n----------------------------------------------------------------\n---- The new charge is being set to a non-standard value -------\n---- This is essentially a fractionally charged system!! -------\n----------------------------------------------------------------\n");

            // change values
            unit_charge = new_charge;
            a_B = (4.0 * Math.PI * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar) / (mass * unit_charge * unit_charge);
            Ry = (mass * unit_charge * unit_charge * unit_charge * unit_charge) / (32.0 * Math.PI * Math.PI * Physics_Base.epsilon * Physics_Base.epsilon * Physics_Base.hbar * Physics_Base.hbar);
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
                return Perform_DoS_Deriv_Integral(band_edge, no_kb_T, alpha, dos_integrand);
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
                return Perform_DoS_Deriv_Integral(band_edge, no_kb_T, alpha, dos_integrand);
            }
        }

        double Perform_DoS_Integral(double band_edge, double no_kb_T, double alpha, OneVariableFunction dos_integrand)
        {
            dos_integrand.Integrator = new GaussKronrodIntegrator();

            // check to see if the band edge is less than a large tolerance... it should be
            // physically, it is difficult for a band edge to dip far below the fermi surface as the charge flowing into the system should pin it
            // here, we use a value of 0.1eV... integral_tol should never be greater than 1.0eV as this is roughly the band gap and holes/electrons
            // should be induced
            double integral_lower = -100.0;
            if (band_edge > integral_lower)
                integral_lower = band_edge;

            double result = alpha * dos_integrand.Integrate(integral_lower, no_kb_T * Physics_Base.kB * temperature);
            return result;
        }

        double Perform_DoS_Deriv_Integral(double band_edge, double no_kb_T, double alpha, OneVariableFunction dos_integrand)
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

        public double Unit_Charge
        { get { return unit_charge; } }
        public double Mass
        { get { return mass; } }

        public void Close()
        {
            Console.WriteLine("Closing density solver");
        }

        #region dftpot
        protected Band_Data dft_pot;
        protected double alpha_dft = 0.1;
        public void Print_DFT_diff(SpinResolved_Data car_dens)
        {
            Console.WriteLine("Maximum absolute difference in DFT potentials = " +  DFT_diff(car_dens).InfinityNorm().ToString("F6"));
        }
        public Band_Data DFT_diff(SpinResolved_Data car_dens)
        {
            if (alpha_dft == 0.0)
                return 0.0 * dft_pot;
            else
                return dft_pot - Get_XC_Potential(car_dens);
        }
        public void Update_DFT_Potential(SpinResolved_Data car_dens)
        {
            if (this.dft_pot == null && alpha_dft != 0.0)
                this.dft_pot = Get_XC_Potential(car_dens);
            else if (alpha_dft == 0.0)
                this.dft_pot = 0.0 * car_dens.Spin_Summed_Data.DeepenThisCopy();
            else
                this.dft_pot = (1.0 - alpha_dft) * dft_pot + alpha_dft * Get_XC_Potential(car_dens);
        }
        public void Reset_DFT_Potential()
        {
            this.dft_pot = null;
        }
        public Band_Data DFT_Potential
        {
            get { return dft_pot; }
        }
        public double DFT_Mixing_Parameter
        {
            get { return alpha_dft; }
            set { alpha_dft = value; }
        }
        #endregion

        #region dft_functions
        double Get_Rs(double charge_density)
        {
            if (charge_density == 0.0)
                return double.PositiveInfinity;

            charge_density = Math.Abs(charge_density / Physics_Base.q_e);    // DFT works on density, not charge 
            return Math.Pow(0.75 / (Math.PI * charge_density), 1.0 / 3.0) / a_B;    // and convert to dimensionless constant
        }
        /*
                static Dictionary<string, double> Unpolarized_Parameters = new Dictionary<string, double>
                {
                    {"kappa", 0.9164},
                    {"A", 0.0622},
                    {"B", -0.096},
                    {"gamma", -0.2846},
                    {"beta1", 1.0529},
                    {"beta2", 0.3334},
                    {"C", 0.0040},
                    {"D", -0.0232}
                };
                static Dictionary<string, double> Polarized_Parameters = new Dictionary<string, double>
                {
                    {"kappa", 1.1546},
                    {"A", 0.03110},
                    {"B", -0.0538},
                    {"gamma", -0.1686},
                    {"beta1", 1.3981},
                    {"beta2", 0.2611},
                    {"C", 0.0014},
                    {"D", -0.0096}
                };*/
        Dictionary<string, double> Unpolarized_Parameters = new Dictionary<string, double>
        {
            {"A", 0.031091},
            {"B", 0.21370},
            {"C", 7.5957},
            {"D", 3.5876},
            {"E", 1.6382},
            {"F", 0.49294}
        };
        Dictionary<string, double> Polarized_Parameters = new Dictionary<string, double>
        {
            {"A", 0.015545},
            {"B", 0.20548},
            {"C", 14.1189},
            {"D", 6.1977},
            {"E", 3.3662},
            {"F", 0.62517}
        };
        Dictionary<string, double> Minus_Spin_Stiffness_Parameters = new Dictionary<string, double>
        {
            {"A", 0.016887},
            {"B", 0.11125},
            {"C", 10.357},
            {"D", 3.6231},
            {"E", 0.88026},
            {"F", 0.49671}
        };

        /// <summary>
        /// returns Perdew-Wang generalised function (see Guiliani & Vignale pp.49-50)
        /// </summary>
        double Get_PW_G_Function(double r_s, Dictionary<string, double> p)
        {
            if (r_s == 0.0)
                return 0.0;
            else if (r_s < 0.0)
                throw new ArgumentOutOfRangeException("Error - r_s must be positive");
            else
            {
                double g = p["A"] * (1.0 + p["C"] * Math.Pow(r_s, 0.5) + p["D"] * r_s + p["E"] * Math.Pow(r_s, 1.5) + p["F"] * r_s * r_s);
                return -4.0 * p["A"] * (1.0 + p["B"] * r_s) * Math.Log(1.0 + 0.5 / g);
            }
        }
        double Get_PW_G_Function_Deriv(double r_s, Dictionary<string, double> p)
        {
            if (r_s == 0.0)
                return 0.0;
            else if (r_s < 0.0)
                throw new ArgumentOutOfRangeException("Error - r_s must be positive");
            else
            {
                double g = p["A"] * (1.0 + p["C"] * Math.Pow(r_s, 0.5) + p["D"] * r_s + p["E"] * Math.Pow(r_s, 1.5) + p["F"] * r_s * r_s);
                double dg = p["A"] * (0.5 * p["C"] * Math.Pow(r_s, -0.5) + p["D"] + 1.5 * p["E"] * Math.Pow(r_s, 0.5) + 2.0 * p["F"] * r_s);

                return -4.0 * p["A"] * p["B"] * Math.Log(1.0 + 0.5 / g)
                        - 4.0 * p["A"] * (1.0 + p["B"] * r_s) / (1.0 + 0.5 / g) * -0.5 * dg / (g * g);
            }
        }

        /// <summary>
        /// return the Perdew - Wang correlation potential
        /// </summary>
        double Get_LSDA_Ec_Potential(double dens_up, double dens_down)
        {
            double polarisation = (dens_up - dens_down) / (dens_up + dens_down);

            double ec, d_ec;                                                 // correlation potential
            double r_s = Get_Rs(dens_up + dens_down);

            if (double.IsInfinity(r_s))
                return 0.0;

            // correlation energy per particle from Guiliani & Vignale (Perdew-Wang form)
            ec = Get_PW_G_Function(r_s, Unpolarized_Parameters);
            if (polarisation != 0.0)
                ec += -1.0 * Get_PW_G_Function(r_s, Minus_Spin_Stiffness_Parameters) * Polarisation_Interpolation(polarisation) / 1.709921 * (1.0 - Math.Pow(polarisation, 4.0))
                        + (Get_PW_G_Function(r_s, Polarized_Parameters) - Get_PW_G_Function(r_s, Unpolarized_Parameters)) * Polarisation_Interpolation(polarisation) * Math.Pow(polarisation, 4.0);

            // and its derivative
            d_ec = Get_PW_G_Function_Deriv(r_s, Unpolarized_Parameters);
            if (polarisation != 0.0)
                d_ec += -1.0 * Get_PW_G_Function_Deriv(r_s, Minus_Spin_Stiffness_Parameters) * Polarisation_Interpolation(polarisation) / 1.709921 * (1.0 - Math.Pow(polarisation, 4.0))
                                + (Get_PW_G_Function_Deriv(r_s, Polarized_Parameters) - Get_PW_G_Function_Deriv(r_s, Unpolarized_Parameters)) * Polarisation_Interpolation(polarisation) * Math.Pow(polarisation, 4.0);

            // potential from LDA form
            double v_c = ec - r_s * d_ec / 3.0;

            // convert from Ry to meV
            v_c *= Ry;

            return v_c;
        }

        double Polarisation_Interpolation(double p)
        {
            return (Math.Pow(1.0 - p, 4.0 / 3.0) + Math.Pow(1.0 + p, 4.0 / 3.0) - 2.0) / (Math.Pow(2.0, 4.0 / 3.0) - 2.0);
        }

        protected double Get_LSDA_Ex_Potential(double dens_up, double dens_down)
        {
            double polarisation = (dens_up - dens_down) / (dens_up + dens_down);
            return 0.5 * Get_Ex_Potential(dens_up + dens_down, Unpolarized_Parameters) * (Math.Pow(1.0 + polarisation, 4.0 / 3.0) + Math.Pow(1.0 - polarisation, 4.0 / 3.0));
        }

        /// <summary>
        /// return the exchange potential (see Guiliani & Vignale p. 34)
        /// </summary>
        double Get_Ex_Potential(double charge_density, Dictionary<string, double> xc_params)
        {
            double v_x;                                                 // exchange potential
            double r_s = Get_Rs(charge_density);

            if (double.IsInfinity(r_s))
                return 0.0;

            // exchange energy per particle for uniform electron gas from Guiliani & Vignale p. 34
            // e_x = - 0.9164 / r_s
            v_x = -1.0 * (4.0 * 0.9164 / 3.0) / r_s;

            // convert from Ry to meV
            v_x *= Ry;

            return v_x;
        }



        /// <summary>
        /// return the derivative, with respect to the density, of the Perdew - Zunger exchange correlation potential
        /// </summary>
        protected double Get_XC_Potential_Deriv(double charge_density)
        {
            double d_vc, d_vx;                                                 // correlation and exchange potentials
            double r_s = Get_Rs(charge_density);

            double density = Math.Abs(charge_density / Physics_Base.q_e);
            double d_rs = -1.0 * r_s / (3.0 * density); // derivative of rs with respect to the density

            if (double.IsInfinity(r_s))
                return 0.0;

            // exchange energy per particle for uniform electron gas from Parr and Yang (Density Functional Theory of Atoms and Molecules)
            // e_x = - 0.9164 / r_s
            d_vx = (4.0 * 0.9164 / 3.0) / (r_s * r_s);

            // correlation energy per particle from Perdew and Zunger (1981) App C, v_xc = (1 - r_s/3 d/dr_s) e_xc (in Ry, note that PZ-1981 is in Ha)
            // e_c = -0.2846 / (1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s)                           r_s > 1
            // e_c = 0.0622 * Math.Log(r_s) - 0.0960 - 0.0232 * r_s + 0.0040 * r_s * Math.Log(r_s)      r_s < 1

            if (r_s > 1)
                d_vc = -0.2846 * (1 + 7.0 * 1.0529 / 6.0 * (0.5 / Math.Sqrt(r_s)) + 4.0 * 0.3334 / 3.0) * Math.Pow(1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s, -2.0)
                                + -0.2846 * (1 + 7.0 * 1.0529 / 6.0 * Math.Sqrt(r_s) + 4.0 * 0.3334 / 3.0 * r_s) * -2.0 * (1.0529 * 0.5 / Math.Sqrt(r_s) + 0.3334) * Math.Pow(1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s, -3.0);
            else
                d_vc = 0.0622 / r_s + (2.0 * 0.0040 / 3.0) * (1.0 + Math.Log(r_s)) + (2.0 * -0.0232 - 0.0040) / 3.0;

            // convert from Ry to meV
            d_vx *= Ry;
            d_vc *= Ry;

            return d_rs * (d_vx + d_vc);
        }

        protected double Get_LSDA_XC_Potential(double dens_up, double dens_down)
        {
            if (dens_up + dens_down == 0.0)
                return 0.0;

            return Get_LSDA_Ex_Potential(dens_up, dens_down) + Get_LSDA_Ec_Potential(dens_up, dens_down);
        }

        public double Get_XC_Potential(double charge_density)
        {
            if (charge_density == 0.0)
                return 0.0;

            // by default, give the unpolarized exchange-correlation potential
            return Get_LSDA_Ex_Potential(0.5 * charge_density, 0.5 * charge_density) + Get_LSDA_Ec_Potential(0.5 * charge_density, 0.5 * charge_density);
        }

        public Band_Data Get_XC_Potential(SpinResolved_Data charge_density)
        {
            Band_Data result;
            Band_Data charge_dens_spin_summed = charge_density.Spin_Summed_Data;
            int dim = charge_dens_spin_summed.Dimension;

            if (dim == 1)
            {
                int nx = charge_dens_spin_summed.vec.Length;
                result = new Band_Data(nx, 0.0);
                for (int i = 0; i < nx; i++)
                    result.vec[i] = Get_XC_Potential(charge_dens_spin_summed.vec[i]);

                return result;
            }
            else if (dim == 2)
            {
                int nx = charge_dens_spin_summed.mat.Rows;
                int ny = charge_dens_spin_summed.mat.Cols;
                result = new Band_Data(nx, ny, 0.0);
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result.mat[i, j] = Get_XC_Potential(charge_dens_spin_summed.mat[i, j]);

                return result;
            }
            else if (dim == 3)
            {
                int nx = charge_dens_spin_summed.vol[0].Rows;
                int ny = charge_dens_spin_summed.vol[0].Cols;
                int nz = charge_dens_spin_summed.vol.Length;
                result = new Band_Data(nx, ny, nz, 0.0);
                for (int k = 0; k < nz; k++)
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                            result.vol[k][i, j] = Get_XC_Potential(charge_dens_spin_summed.vol[k][i, j]);

                return result;
            }
            else
                throw new NotImplementedException();
        }

        public Band_Data Get_XC_Potential_Deriv(SpinResolved_Data charge_density)
        {
            Band_Data result;
            Band_Data charge_dens_spin_summed = charge_density.Spin_Summed_Data;
            int dim = charge_dens_spin_summed.Dimension;

            if (dim == 1)
            {
                int nx = charge_dens_spin_summed.vec.Length;
                result = new Band_Data(new DoubleVector(nx));
                for (int i = 0; i < nx; i++)
                    result.vec[i] = Get_XC_Potential_Deriv(charge_dens_spin_summed.vec[i]);

                return result;
            }
            else if (dim == 2)
            {
                int nx = charge_dens_spin_summed.mat.Rows;
                int ny = charge_dens_spin_summed.mat.Cols;
                result = new Band_Data(new DoubleMatrix(nx, ny));
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result.mat[i, j] = Get_XC_Potential_Deriv(charge_dens_spin_summed.mat[i, j]);

                return result;
            }
            else if (dim == 3)
            {
                throw new NotImplementedException();
            }
            else
                throw new NotImplementedException();
        }
        #endregion
    }
}
