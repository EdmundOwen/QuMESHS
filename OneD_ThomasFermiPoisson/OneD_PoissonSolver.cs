using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace OneD_ThomasFermiPoisson
{
    class OneD_PoissonSolver : Potential_Base
    {
        // parameters for regular grid solve
        int nz;
        double dz;

        double top_bc, bottom_bc;

        // 
        public OneD_PoissonSolver(double dz, int nz, double top_bc, double bottom_bc)
        {
            this.nz = nz; this.dz = dz;

            this.top_bc = top_bc; this.bottom_bc = bottom_bc;
        }

        public DoubleVector Get_Potential(DoubleVector density)
        {
            // calculate potential on a regular grid (not ideal, or scalable)
            return Get_Potential_On_Regular_Grid(density);
        }

        DoubleVector Get_Potential_On_Regular_Grid(DoubleVector density)
        {
            // generate Laplacian matrix (spin-resolved)
            DoubleMatrix laplacian = Generate_Laplacian();

            // set boundary conditions
            density[0] = top_bc; density[nz] = top_bc;
            density[nz - 1] = bottom_bc; density[2 * nz - 1] = bottom_bc;

            DoubleVector potential = NMathFunctions.Solve(laplacian, density);

            return potential;
        }

        /// <summary>
        /// Generates a spin-resolved laplacian matrix in one-dimension on a regular grid with Dirichlet BCs
        /// </summary>
        /// <returns></returns>
        DoubleMatrix Generate_Laplacian()
        {
            // the factor which multiplies the Laplace equation
            double factor = -1.0 * epsilon / (dz * dz);

            DoubleMatrix result = new DoubleMatrix(2 * nz, 2 * nz);
            for (int i = 0; i < 2 * nz - 1; i++)
            {
                // on-diagonal term
                result[i, i] = 2.0 * factor;
                // off-diagonal
                result[i + 1, i] = -1.0 * factor;
                result[i, i + 1] = -1.0 * factor;
            }
            // and delete spin-sector coupling terms 
            result[nz - 1, nz] = 0.0;
            result[nz, nz - 1] = 0.0;

            // and fix boundary conditions for spin-up
            result[0, 0] = 1.0 * factor;
            result[0, 1] = 0.0;
            result[nz - 1, nz - 1] = 1.0 * factor;
            result[nz - 1, nz - 2] = 0.0;

            // and similarly for spin-down
            result[nz, nz] = 1.0 * factor;
            result[nz, nz + 1] = 0.0;
            result[2 * nz - 1, 2 * nz - 1] = 1.0 * factor;
            result[2 * nz - 1, 2 * nz - 2] = 0.0;

            return result;
        }
    }
}
