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
        public const double mass = 0.067 * 5.68562958e-3;       // (meV) (ps)^2 (nm)^-2 with GaAs effective mass
        public const double epsilon_0 = 1.41859713;             // (zC)^2 (nm)^-1 (meV)^-1 for vacuum
        public const double epsilon_r = 13.0;                   // relative permittivity for vacuum -> GaAs
        public const double epsilon = epsilon_0 * epsilon_r;    // (zC)^2 (nm)^-1 (meV)^-1 for GaAs
        public const double kB = 0.086173324;                   // (meV) (K)^-1
        
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
        /// gets 3D spin-resolved density of states for given potential and energy
        /// </summary>
        public static double Get_3D_DensityofStates(double energy, double potential)
        {
            // if the energy is below the potential (i.e. below the conduction band) return zero
            if (energy < potential)
                return 0.0;

            // generate density of states
            return Math.Pow(2.0 * mass * mass * mass * (energy - potential), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);

            // calculate dk^2 / dE and k^2
            //double dk2_dE = (2 * mass) / hbar * hbar;
            //double k2 = (energy - potential) * dk2_dE;

            // prefactor for number of states per unit volume... ie (4 pi / 3) / (2 pi)^3
            //double geometric_prefactor = 1.0 / (6 * Math.PI * Math.PI);

            //double density_of_states = geometric_prefactor * dk2_dE * 1.5 * Math.Pow(k2, 0.5);
            //return density_of_states;
        }
    }
}
