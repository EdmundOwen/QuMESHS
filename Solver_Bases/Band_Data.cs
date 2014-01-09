using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    /// <summary>
    /// a class for holding band structure data which can be either
    /// 1, 2, or 3 dimensional
    /// </summary>
    public class Band_Data
    {
        int dim;
        public DoubleVector vec;
        public DoubleMatrix mat;
        public DoubleMatrix[] vol;

        public Band_Data(DoubleVector vec)
        {
            this.vec = vec;
            dim = 1;
        }

        public Band_Data(DoubleMatrix mat)
        {
            this.mat = mat;
            dim = 2;
        }

        public Band_Data(DoubleMatrix[] vol)
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

            set
            {
                if (Dimension == 1)
                    vec[i] = value;
                else if (Dimension == 2)
                {
                    int x = i % mat.Cols;
                    int y = (int)((i - x) / mat.Cols);
                    mat[y, x] = value;
                }
                else if (Dimension == 3)
                {
                    int matsize = vol[0].Rows * vol[0].Cols;
                    int index1 = i % matsize;
                    int x = i % vol[0].Cols;
                    int y = (int)((i - x) / vol[0].Cols);
                    int z = (int)((i - index1) / matsize);

                    vol[z][y, x] = value;
                }
                else
                    throw new NotImplementedException();
            }
        }

        public static Band_Data operator +(Band_Data vec1, Band_Data vec2)
        {
            if (vec1.Dimension != vec2.Dimension)
                throw new RankException();

            if (vec1.Dimension == 1)
                return new Band_Data(vec1.vec + vec2.vec);
            else if (vec1.Dimension == 2)
                return new Band_Data(vec1.mat + vec2.mat);
            else if (vec1.Dimension == 3)
            {
                Band_Data result = new Band_Data(new DoubleMatrix[vec1.vol.Length]);
                for (int i = 0; i < vec1.vol.Length; i++)
                    result.vol[i] = vec1.vol[i] + vec2.vol[i];

                return result;
            }
            else
                throw new NotImplementedException();
        }

        public static Band_Data operator -(Band_Data vec1, Band_Data vec2)
        {
            if (vec1.Dimension != vec2.Dimension)
                throw new RankException();

            if (vec1.Dimension == 1)
                return new Band_Data(vec1.vec - vec2.vec);
            else if (vec1.Dimension == 2)
                return new Band_Data(vec1.mat - vec2.mat);
            else if (vec1.Dimension == 3)
            {
                Band_Data result = new Band_Data(new DoubleMatrix[vec1.vol.Length]);
                for (int i = 0; i < vec1.vol.Length; i++)
                    result.vol[i] = vec1.vol[i] - vec2.vol[i];

                return result;
            }
            else
                throw new NotImplementedException();
        }

        public static Band_Data operator *(double scalar, Band_Data data)
        {
            for (int i = 0; i < data.Length; i++)
                data[i] *= scalar;

            return data;
        }

        public static Band_Data operator *(Band_Data data, double scalar)
        {
            return scalar * data;
        }

        public static Band_Data operator /(Band_Data data, double scalar)
        {
            return (1.0 / scalar) * data;
        }

        public void Save_1D_Data(string filename, double dz, double zmin)
        {
            // check that the dimension of the density is correct
            if (this.Dimension != 1)
                throw new RankException();

            // open stream
            StreamWriter sw = new StreamWriter(filename);

            int nz = this.vec.Length;

            // save out positions
            sw.WriteLine("x " + nz.ToString());
            for (int i = 0; i < nz; i++)
                sw.WriteLine(((float)(i * dz + zmin)).ToString());

            // save out densities
            sw.WriteLine();
            sw.WriteLine("data");
            for (int i = 0; i < nz; i++)
                sw.WriteLine(((float)this.vec[i]).ToString());

            sw.Close();
        }
            
        public void Save_2D_Data(string filename, double dy, double dz, double ymin, double zmin)
        {
            // check that the dimension of the density is correct
            if (this.Dimension != 2)
                throw new RankException();

            // open stream
            StreamWriter sw = new StreamWriter(filename);

            int ny = this.mat.Rows;
            int nz = this.mat.Cols;

            // save out positions
            sw.WriteLine("x " + ny.ToString());
            for (int i = 0; i < ny; i++)
                //for (int j = 0; j < nz; j++)
                sw.Write((i * dy + ymin).ToString() + '\t');

            sw.WriteLine();
            sw.WriteLine();
            sw.WriteLine("y " + nz.ToString());
            for (int j = 0; j < nz; j++)
                sw.Write(((float)(j * dz + zmin)).ToString() + '\t');

            // save out densities
            sw.WriteLine();
            sw.WriteLine("data");
            for (int i = 0; i < nz; i++)
            {
                for (int j = 0; j < ny; j++)
                    if (Math.Abs(this.mat[j, i]) < 1e-20)
                        sw.Write("0\t");
                    else
                        // note that the ordering is y first, then z -- this is FlexPDE specific
                        sw.Write(((float)this.mat[j, i]).ToString() + '\t');
                sw.WriteLine();
            }

            sw.Close();
        }

        public void Save_3D_Data(string filename, double dx, double dy, double dz, double xmin, double ymin, double zmin)
        {
            throw new NotImplementedException();
        }

        public int Dimension
        {
            get { return dim; }
        }

        /// <summary>
        /// returns the number of data points in the Band_Data object
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
