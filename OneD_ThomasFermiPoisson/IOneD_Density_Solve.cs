using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;

namespace OneD_ThomasFermiPoisson
{
    interface IOneD_Density_Solve : IDensity_Solve
    {
        double Zmin_Pot { set; }
        double Dz_Pot { set; }
    }
}
