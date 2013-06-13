using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public abstract class Density_Solver
    {
        protected double fermi_Energy, temperature;
        protected int nx, ny, nz;

        // physical constants
        protected const double hbar = 0.658211814;             // (meV) (ps)
        protected const double q_e = -0.160217646;             // (aC)
        protected const double mass = 0.067 * 5.68562958e-3;   // (meV) (ps)^2 (nm)^-2 with GaAs effective mass
        protected const double epsilon = 13 * 1.41859713e-6;   // (aC)^2 (nm)^-1 (meV)^-1 for GaAs
        protected const double kB = 0.086173324;               // (meV) (K)^-1

        public Density_Solver(double fermi_Energy, double temperature, int nx, int ny, int nz)
        {
            this.fermi_Energy = fermi_Energy; this.temperature = temperature;
            this.nx = nx; this.ny = ny; this.nz = nz;
        }

        /// <summary>
        /// Gets the dopent occupation probability at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Dopent_Occupation(double energy)
        {
            return Get_Dopent_Occupation(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Calculates the dopent occupation probability for arbitrary energy, E_f and T
        /// </summary>
        protected double Get_Dopent_Occupation(double energy, double E_f, double T)
        {
            if (T == 0)
                if (energy > E_f)
                    return 0.0;
                else
                    return 1.0;
            else
                return 2.0 / (Math.Exp((energy - fermi_Energy) / (kB * T)) + 2.0);
        }

        /// <summary>
        /// Gets the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function(double energy)
        {
            return Get_Fermi_Function(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Calculates the fermi function for arbitrary energy, E_f and T
        /// </summary>
        protected double Get_Fermi_Function(double energy, double E_f, double T)
        {
            if (T == 0)
                if (energy > E_f)
                    return 0.0;
                else
                    return 1.0;
            else
                return 1.0 / (Math.Exp((energy - fermi_Energy) / (kB * T)) + 1.0);
        }
    }
}
