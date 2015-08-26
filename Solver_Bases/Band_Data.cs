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
        string value_location = "";

        Band_Data laplacian_band_data;
        

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

        public Band_Data(int nz, double val)
        { 
            this.vec = new DoubleVector(nz, val);
            dim = 1;
        }

        public Band_Data(int ny, int nz, double val)
        {
            this.mat = new DoubleMatrix(ny, nz, val);
            dim = 2;
        }

        public Band_Data(int nx, int ny, int nz, double val)
        {
            this.vol = new DoubleMatrix[nz];
            for (int i = 0; i < nz; i++)
                this.vol[i] = new DoubleMatrix(nx, ny, val);
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
                    int x = index1 % vol[0].Cols;
                    int y = (int)((index1 - x) / vol[0].Cols);
                    int z = (int)((i - index1) / matsize);

                    return (double)vol[z][y, x];
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
                    int y = (int)((i - x) / vol[0].Cols) % vol[0].Rows;
                    int z = (int)((i - index1) / matsize);

                    vol[z][y, x] = value;
                }
                else
                    throw new NotImplementedException();
            }
        }

        public static Band_Data operator +(Band_Data vec1, Band_Data vec2)
        {
            Band_Data result;

            // return null if both of the inputs are null
            if (vec1 == null || vec2 == null)
                return null;

            if (vec1.Dimension != vec2.Dimension)
                throw new RankException();

            if (vec1.Dimension == 1)
                result = new Band_Data(vec1.vec + vec2.vec);
            else if (vec1.Dimension == 2)
                result = new Band_Data(vec1.mat + vec2.mat);
            else if (vec1.Dimension == 3)
            {
                result = new Band_Data(new DoubleMatrix[vec1.vol.Length]);
                for (int i = 0; i < vec1.vol.Length; i++)
                    result.vol[i] = vec1.vol[i] + vec2.vol[i];
            }
            else
                throw new NotImplementedException();

            result.Laplacian = vec1.Laplacian + vec2.Laplacian;
            return result;
        }

        public static Band_Data operator -(Band_Data vec1, Band_Data vec2)
        {
            Band_Data result;

            // return null if both of the inputs are null
            if (vec1 == null || vec2 == null)
                return null;

            if (vec1.Dimension != vec2.Dimension)
                throw new RankException();

            if (vec1.Dimension == 1)
                result = new Band_Data(vec1.vec - vec2.vec);
            else if (vec1.Dimension == 2)
                result = new Band_Data(vec1.mat - vec2.mat);
            else if (vec1.Dimension == 3)
            {
                result = new Band_Data(new DoubleMatrix[vec1.vol.Length]);
                for (int i = 0; i < vec1.vol.Length; i++)
                    result.vol[i] = vec1.vol[i] - vec2.vol[i];
            }
            else
                throw new NotImplementedException();

            result.Laplacian = vec1.Laplacian - vec2.Laplacian;
            return result;
        }

        public static Band_Data operator *(double scalar, Band_Data data)
        {
            Band_Data result;

            // return null if the inputs is null
            if (data == null)
                return null;

            if (data.Dimension == 1)
                result = new Band_Data(scalar * data.vec);
            else if (data.Dimension == 2)
                result = new Band_Data(scalar * data.mat);
            else if (data.Dimension == 3)
            {
                result = new Band_Data(new DoubleMatrix[data.vol.Length]);
                for (int i = 0; i < data.vol.Length; i++)
                    result.vol[i] = scalar * data.vol[i];
            }
            else
                throw new NotImplementedException();

            result.Laplacian = scalar * data.Laplacian;
            return result;
        }

        public static Band_Data operator *(Band_Data data, double scalar)
        {
            return scalar * data;
        }

        public static Band_Data operator /(Band_Data data, double scalar)
        {
            return (1.0 / scalar) * data;
        }

        /// <summary>
        /// Initiates a laplacian for this Band_Data class of the correct size
        /// </summary>
        public void Initiate_Laplacian()
        {
            this.laplacian_band_data = 0.0 * this;
        }

        public void Save_1D_Data(string filename, double dz, double zmin)
        {
            this.value_location = filename;
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
            this.value_location = filename;
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
            this.value_location = filename;
            // check that the dimension of the density is correct
            if (this.Dimension != 3)
                throw new RankException();

            // open stream
            StreamWriter sw = new StreamWriter(filename);

            int nz = this.vol.Length;
            int nx = this.vol[0].Rows;
            int ny = this.vol[0].Cols;

            // save out positions
            sw.WriteLine("x " + nx.ToString());
            for (int i = 0; i < nx; i++)
                sw.Write((i * dx + xmin).ToString() + '\t');

            sw.WriteLine();
            sw.WriteLine();
            sw.WriteLine("y " + ny.ToString());
            for (int j = 0; j < ny; j++)
                sw.Write(((float)(j * dy + ymin)).ToString() + '\t');

            sw.WriteLine();
            sw.WriteLine();
            sw.WriteLine("z " + nz.ToString());
            for (int k = 0; k < nz; k++)
                sw.Write(((float)(k * dz + zmin)).ToString() + '\t');

            // save out densities
            sw.WriteLine();
            sw.WriteLine("data");
            for (int i = 0; i < nz; i++)
                for (int j = 0; j < ny; j++)
                {
                    for (int k = 0; k < nx; k++)
                        if (Math.Abs(this.vol[i][k, j]) < 1e-20)
                            sw.Write("0\t");
                        else
                            // note that the ordering is x first, then y, then z -- this is FlexPDE specific
                            sw.Write(((float)this.vol[i][k, j]).ToString() + '\t');
                    sw.WriteLine();
                }

            sw.Close();
        }

        /// <summary>
        /// saves the data contained within the class as an unstructured list of data points
        /// </summary>
        /// <param name="filename"></param>
        public void Save_Data(string filename)
        {
            this.value_location = filename;
            StreamWriter sw = new StreamWriter(filename);
            if (dim == 1)
                for (int i = 0; i < Length; i++)
                    sw.WriteLine(this[i]);
            else if (dim == 2)
                for (int i = 0; i < mat.Cols; i++)
                    for (int j = 0; j < mat.Rows; j++)
                        sw.WriteLine(mat[j, i]);
            else if (dim == 3)
                for (int k = 0; k < vol.Length; k++ )
                    for (int i = 0; i < vol[0].Cols; i++)
                        for (int j = 0; j < vol[0].Rows; j++)
                            sw.WriteLine(vol[k][j, i]);
            else
                throw new NotImplementedException();
            sw.Close();
        }

        public static Band_Data Parse_Band_Data(string location, string[] input_data, int nz)
        {
            // and check that there is the right number of data points back
            if (input_data.Length != nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleVector(nz));
            for (int i = 0; i < nz; i++)
                result.vec[i] = double.Parse(input_data[i]);

            result.value_location = location;
            return result;
        }

        public static Band_Data Parse_Band_Data(string location, string[] input_data, int ny, int nz)
        {
            // and check that there is the right number of data points back
            if (input_data.Length != ny * nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleMatrix(ny, nz));
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nz; j++)
                    result.mat[i, j] = double.Parse(input_data[j * ny + i]);
            }

            result.value_location = location;
            return result;
        }

        public static Band_Data Parse_Band_Data(string location, string[] input_data, int nx, int ny, int nz)
        {
            // and check that there is the right number of data points back
            if (input_data.Length != nx * ny * nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(nx, ny, nz, 0.0);
            for (int k = 0; k < nz; k++)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result.vol[k][i, j] = double.Parse(input_data[k * nx * ny + j * nx + i]);

            result.value_location = location;
            return result;
        }

        public Band_Data DeepenThisCopy()
        {
            if (dim == 1)
            {
                DoubleVector result = new DoubleVector(Length);
                for (int i = 0; i < Length; i++)
                    result[i] = this[i];

                return new Band_Data(result);
            }
            else if (dim == 2)
            {
                DoubleMatrix result = new DoubleMatrix(mat.Rows, mat.Cols);
                for (int i = 0; i < mat.Rows; i++)
                    for (int j = 0; j < mat.Cols; j++)
                        result[i, j] = mat[i, j];

                return new Band_Data(result);
            }
            else if (dim == 3)
            {
                int nx = vol[0].Rows;
                int ny = vol[0].Cols;
                int nz = vol.Length;

                Band_Data result = new Band_Data(nx, ny, nz, 0.0);

                for (int k = 0; k < nz; k++)
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                            result.vol[k][i, j] = this.vol[k][i, j];

                return result;
            }
            else
                throw new NotImplementedException();
        }

        public double Max()
        {
            if (dim == 1)
                return vec.Max();
            else if (dim == 2)
                return mat.Max();
            else if (dim == 3)
            {
                DoubleVector mat_max = new DoubleVector(vol.Length);
                for (int i = 0; i < vol.Length; i++)
                    mat_max[i] = vol[i].Max();

                return mat_max.Max();
            }
            else
                throw new NotImplementedException();
        }

        public double Min()
        {
            if (dim == 1)
                return vec.Min();
            else if (dim == 2)
                return mat.Min();
            else if (dim == 3)
            {
                DoubleVector mat_min = new DoubleVector(vol.Length);
                for (int i = 0; i < vol.Length; i++)
                    mat_min[i] = vol[i].Min();

                return mat_min.Min();
            }
            else
                throw new NotImplementedException();
        }

        public double InfinityNorm()
        {
            if (dim == 1)
                return vec.InfinityNorm();
            else if (dim == 2)
            {
                DoubleVector mat_absmax = new DoubleVector(mat.Rows);
                for (int i = 0; i < mat.Rows; i++)
                    mat_absmax[i] = mat.Row(i).InfinityNorm();

                return mat_absmax.InfinityNorm();
            }
            else if (dim == 3)
            {
                DoubleVector vol_absmax = new DoubleVector(vol.Length);
                for (int i = 0; i < vol.Length; i++)
                    vol_absmax[i] = vol[i].InfinityNorm();

                return vol_absmax.InfinityNorm();
            }
            else
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

        public string File_Location
        {
            get { return value_location; }
        }

        public Band_Data Laplacian
        {
            get 
            {
                //if (laplacian_band_data == null)
                    //throw new TypeInitializationException("laplacian_band_data", new Exception("Error - Band data class for the laplacian has not been calculated"));
                return laplacian_band_data; 
            }
            set { laplacian_band_data = value; }
        }
    }
}
