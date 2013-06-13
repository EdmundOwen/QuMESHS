using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public abstract class Potential_Base
    {
        // physical constants
        protected const double q_e = -0.160217646;             // (aC)
        protected const double epsilon = 13 * 1.41859713e-6;   // (aC)^2 (nm)^-1 (meV)^-1 for GaAs
    }
}
