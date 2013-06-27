using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public class Input_Band_Structure
    {
        public static DoubleVector GetBandStructure(string filename)
        {
            // temporary band structure... (spin-degenerate)
            DoubleVector result = new DoubleVector(nz, 2100.0);
            for (int k = 100; k < 130; k++)
            {
                //int k = 10;
                result[k] = 1400.0;
            }

            return result;
        }

        /// <summary>
        /// returns a DoubleMatrix with the given band structure planarised in the transverse direction
        /// </summary>
        public static DoubleMatrix Expand_BandStructure(DoubleVector structure, int ny)
        {
            DoubleMatrix result = new DoubleMatrix(ny, structure.Length);
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < structure.Length; j++)
                    result[i, j] = structure[j];

            return result;
        }

        /// <summary>
        /// returns a DoubleMatrix with the given band structure planarised in the transverse direction
        /// </summary>
        public static DoubleMatrix[] Expand_BandStructure(DoubleVector structure, int nx, int ny)
        {
            DoubleMatrix[] result = new DoubleMatrix[nx];

            for (int i = 0; i < nx; i++)
            {
                result[i] = new DoubleMatrix(ny, structure.Length);
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < structure.Length; k++)
                        result[i][j, k] = structure[k];
            }

            return result;
        }
    }
}
