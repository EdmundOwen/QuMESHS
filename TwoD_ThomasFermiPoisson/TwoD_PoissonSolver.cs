using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_PoissonSolver
    {
        public TwoD_PoissonSolver()
        {
        }

        public DoubleVector Get_Well_Potential(int nx)
        {
            Console.WriteLine("THIS IS A PROTOTYPE! NO POISSON SOLVER IMPLEMENTED!!!");

            DoubleVector result = new DoubleVector(nx);
            for (int i = 0; i < nx; i++)
                result[i] = (i - (double)(nx - 1) / 2) * (i - (double)(nx - 1) / 2);

            return result;
        }
    }
}
