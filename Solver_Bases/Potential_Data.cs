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

        /// <summary>
        /// returns a value for this data type.
        /// Note that value ordering is not very well defined
        /// </summary>
        public double this[int i]
        {
            get 
            {
                if (Dimension == 1)
                    return vec[i];
                else if (Dimension == 2)
                {
                    int x = i % mat.Cols;
                    int y = (int)((i - x) / mat.Cols);
                    return (double)mat[y, x];
                }
                else if (Dimension == 3)
                {
                    int matsize = vol[0].Rows * vol[0].Cols;
                    int index1 = i % matsize;
                    int x = i % vol[0].Cols;
                    int y = (int)((i - x) / vol[0].Cols);
                    int z = (int)((i - index1) / matsize);

                    return (double)vol[z][y,x];
                }
                else
                    throw new NotImplementedException();
            }
        }

        public static Potential_Data operator +(Potential_Data vec1, Potential_Data vec2)
        {
            if (vec1.Dimension != vec2.Dimension)
                throw new RankException();

            if (vec1.Dimension == 1)
                return new Potential_Data(vec1.vec + vec2.vec);
            else if (vec1.Dimension == 2)
                return new Potential_Data(vec1.mat + vec2.mat);
            else if (vec1.Dimension == 3)
            {
                Potential_Data result = new Potential_Data(new DoubleMatrix[vec1.vol.Length]);
                for (int i = 0; i < vec1.vol.Length; i++)
                    result.vol[i] = vec1.vol[i] + vec2.vol[i];

                return result;
            }
            else
                throw new NotImplementedException();
        }

        public static Potential_Data operator -(Potential_Data vec1, Potential_Data vec2)
        {
            if (vec1.Dimension != vec2.Dimension)
                throw new RankException();

            if (vec1.Dimension == 1)
                return new Potential_Data(vec1.vec - vec2.vec);
            else if (vec1.Dimension == 2)
                return new Potential_Data(vec1.mat - vec2.mat);
            else if (vec1.Dimension == 3)
            {
                Potential_Data result = new Potential_Data(new DoubleMatrix[vec1.vol.Length]);
                for (int i = 0; i < vec1.vol.Length; i++)
                    result.vol[i] = vec1.vol[i] - vec2.vol[i];

                return result;
            }
            else
                throw new NotImplementedException();
        }


        public int Dimension
        {
            get { return dim; }
        }

        /// <summary>
        /// returns the number of data points in the Potential_Data object
        /// </summary>
        public int Length
        {
            get
            {
                if (Dimension == 1)
                    return vec.Length;
                else if (Dimension == 2)
                    return mat.Cols * mat.Rows;
                else if (Dimension == 3)
                    return vol.Length * vol[0].Cols * vol[0].Rows;
                else
                    throw new NotImplementedException();
            }
        }
    }
}
