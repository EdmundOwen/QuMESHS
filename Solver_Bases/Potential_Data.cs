using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    /// <summary>
    /// a class for holding data used in the potential which can be either
    /// 1, 2, or 3 dimensional
    /// </summary>
    public class Potential_Data
    {
        int dim;
        public DoubleVector vec;
        public DoubleMatrix mat;
        public DoubleMatrix[] vol;

        public Potential_Data(DoubleVector vec)
        {
            this.vec = vec;
            dim = 1;
        }

        public Potential_Data(DoubleMatrix mat)
        {
            this.mat = mat;
            dim = 2;
        }

        public Potential_Data(DoubleMatrix[] vol)
        {
            this.vol = vol;
            dim = 3;
        }

        public int Dimension
        {
            get { return dim; }
        }
    }
}
