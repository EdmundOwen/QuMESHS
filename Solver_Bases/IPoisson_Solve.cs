using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public interface IPoisson_Solve
    {
        Band_Data Calculate_Laplacian(Band_Data input_vec);
    }
}
