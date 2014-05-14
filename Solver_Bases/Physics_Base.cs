using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public static class Physics_Base
    {
        // physical constants
        public const double hbar = 0.658211814;                 // (meV) (ps)
        public const double q_e = 160.217646;                   // (zC) is positive as in these definitions it is the elementary charge
        public const double m_e = 5.68562958e-3;                // (meV) (ps)^2 (nm)^-2 of the free electron mass
        public const double mass = 0.067 * m_e;                 // (meV) (ps)^2 (nm)^-2 with GaAs effective mass
        public const double energy_V_to_meVpzC = 6.2415093;     // conversion factor from V to meV per zC
        public const double epsilon_0 = 1.41859713;             // (zC)^2 (nm)^-1 (meV)^-1 for vacuum
        public const double epsilon_r = 13.0;                   // relative permittivity for vacuum -> GaAs
        public const double epsilon_r_GaAs = 12.9;              // relative permittivity for vacuum -> GaAs
        public const double epsilon_r_AlGaAs = 12.0;            // relative permittivity for vacuum -> AlGaAs
        public const double epsilon_pmma = 2.6;                 // relative permittivity for vacuum -> PMMA
        public const double epsilon = epsilon_0 * epsilon_r;    // (zC)^2 (nm)^-1 (meV)^-1 for GaAs
        public const double kB = 0.086173324;                   // (meV) (K)^-1
        public const double a_B = (4.0 * Math.PI * epsilon * hbar * hbar) / (mass * q_e * q_e);       // Bohr radius for GaAs in nm
        public const double a_0 = (4.0 * Math.PI * epsilon_0 * hbar * hbar) / (m_e * q_e * q_e);      // True Bohr radius
        public const double alpha = 1.0 / 137.35999074;         // fine structure constant
        public const double Ry_0 = 13605.69253;                 // rydberg in meV in free space
        public const double Ry = (mass * q_e * q_e * q_e * q_e) / (32.0 * Math.PI * Math.PI * epsilon * epsilon * hbar * hbar);     // rydberg in meV in GaAs

        /// <summary>
        /// Calculates the fermi function for arbitrary energy, E_f and T
        /// </summary>
        public static double Get_Fermi_Function(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy > mu)
                    return 0.0;
                else
                    return 1.0;
            else
                return 1.0 / (Math.Exp((energy - mu) / (kB * T)) + 1.0);
        }

        /// <summary>
        /// Calculates the derivative of the fermi function for arbitrary energy, E_f and T
        /// </summary>
        public static double Get_Fermi_Function_Derivative(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy == mu)
                    return 1.0;
                else
                    return 0.0;
            else
            {
                double exponent = (energy - mu) / kB * T;

                if (double.IsInfinity(Math.Exp(exponent)))
                    return 0.0;
                else
                    return Math.Exp(exponent) / (Math.Pow((Math.Exp(exponent) + 1.0), 2.0) * kB * T);
            }
        }

        /// <summary>
        /// Calculates the spin-resolved fermi function for a dopent.
        /// This is different from the typical fermi function as double occupation of the donor is not allowed
        /// </summary>
        public static double Get_Dopent_Fermi_Function(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy > mu)
                    return 0.0;
                else
                    return 1.0;
            else
                return 1.0 / (Math.Exp((energy - mu) / (kB * T)) + 2.0);
        }

        /// <summary>
        /// Calculates the derivative of the spin-resolved fermi function for a dopent.
        /// This is different from the typical fermi function as double occupation of the donor is not allowed
        /// </summary>
        public static double Get_Dopent_Fermi_Function_Derivative(double energy, double mu, double T)
        {
            if (T == 0)
                if (energy == mu)
                    return 1.0;
                else
                    return 0.0;
            else
            {
                double exponent = (energy - mu) / kB * T;

                if (double.IsInfinity(Math.Exp(exponent)))
                    return 0.0;
                else
                    return Math.Exp(exponent) / (Math.Pow((Math.Exp(exponent) + 2.0), 2.0) * kB * T);
            }
        }

        /// <summary>
        /// gets 3D spin-resolved density of states for given potential and energy
        /// </summary>
        public static double Get_Electron_3D_DensityofStates(double energy, double conduction_band_edge)
        {
            // if the energy is below the potential (i.e. below the conduction band) return zero
            if (energy < conduction_band_edge)
                return 0.0;

            // generate density of states
            return Math.Pow(2.0 * mass * mass * mass * (energy - conduction_band_edge), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);

            // calculate dk^2 / dE and k^2
            //double dk2_dE = (2 * mass) / hbar * hbar;
            //double k2 = (energy - potential) * dk2_dE;

            // prefactor for number of states per unit volume... ie (4 pi / 3) / (2 pi)^3
            //double geometric_prefactor = 1.0 / (6 * Math.PI * Math.PI);

            //double density_of_states = geometric_prefactor * dk2_dE * 1.5 * Math.Pow(k2, 0.5);
            //return density_of_states;
        }

        /// <summary>
        /// gets 3D spin-resolved density of states for electrons given potential and energy
        /// </summary>
        public static double Get_Hole_3D_DensityofStates(double energy, double valence_band_edge)
        {
            // if the energy is below the potential (i.e. below the conduction band) return zero
            if (energy > valence_band_edge)
                return 0.0;

            // generate density of states
            return Math.Pow(2.0 * mass * mass * mass * (valence_band_edge - energy), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);
        }

        /// <summary>
        /// return the Perdew - Zunger exchange correlation potential
        /// </summary>
        public static double Get_XC_Potential(double charge_density)
        {
            double v_c, v_x;                                                 // correlation and exchange potentials
            double density = Math.Abs(charge_density) / Physics_Base.q_e;    // DFT works on density, not charge 
            double r_s = Math.Pow(0.75 / (Math.PI * density), 1.0 / 3.0) / Physics_Base.a_B;    // and convert to dimensionless constant

            if (density == 0 || double.IsInfinity(r_s))
                return 0.0;

            // exchange energy per particle for uniform electron gas from Parr and Yang (Density Functional Theory of Atoms and Molecules)
            // e_x = - 0.9164 / r_s
            v_x = -1.0 * (4.0 * 0.9164 / 3.0) / r_s;

            // correlation energy per particle from Perdew and Zunger (1981) App C, v_xc = (1 - r_s/3 d/dr_s) e_xc (in Ry, note that PZ-1981 is in Ha)
            // e_c = -0.2846 / (1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s)                           r_s > 1
            // e_c = 0.0622 * Math.Log(r_s) - 0.0960 - 0.0232 * r_s + 0.0040 * r_s * Math.Log(r_s)      r_s < 1

            if (r_s > 1)
                v_c = -0.2846 * (1 + 7.0 * 1.0529 / 6.0 * Math.Sqrt(r_s) + 4.0 * 0.3334 / 3.0 * r_s) * Math.Pow(1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s, -2.0);
            else
                v_c = 0.0622 * Math.Log(r_s) - (0.0960 - 0.0622 / 3.0) + (2.0 * 0.0040 / 3.0) * r_s * Math.Log(r_s) + (2.0 * -0.0232 - 0.0040) / 3.0 * r_s;

            // convert from Ry to meV
            v_x *= Physics_Base.Ry;
            v_c *= Physics_Base.Ry;

            return v_x + v_c;
        }
    }
}
