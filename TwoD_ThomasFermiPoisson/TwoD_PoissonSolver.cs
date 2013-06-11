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

        public DoubleVector Get_Well_Potential(int ny)
        {
            Console.WriteLine("THIS IS A PROTOTYPE! NO POISSON SOLVER IMPLEMENTED!!!");

            DoubleVector result = new DoubleVector(ny);
            for (int i = 0; i < ny; i++)
                result[i] = (i - (double)(ny - 1) / 2) * (i - (double)(ny - 1) / 2);

            return result;
        }
    }
}
