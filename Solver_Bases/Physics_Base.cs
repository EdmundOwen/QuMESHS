using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public abstract class Physics_Base
    {
        // physical constants
        protected const double hbar = 0.658211814;             // (meV) (ps)
        protected const double q_e = -160.217646;              // (fC)
        protected const double mass = 0.067 * 5.68562958e-3;   // (meV) (ps)^2 (nm)^-2 with GaAs effective mass
        protected const double epsilon = 13 * 1.41859713;      // (fC)^2 (nm)^-1 (meV)^-1 for GaAs
        protected const double kB = 0.086173324;               // (meV) (K)^-1
        
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
                return 1.0 / (Math.Exp((energy - E_f) / (kB * T)) + 1.0);
        }
    }
}
