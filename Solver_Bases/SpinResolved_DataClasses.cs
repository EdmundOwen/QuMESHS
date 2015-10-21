/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{    
    /// <summary>
    /// A class for holding spin resolved DoubleVector classes from CenterSpace types
    /// This is just a pair of DoubleVectors which can be indexed
    /// </summary>
    public class SpinResolved_DoubleVector
    {
        DoubleVector[] spin_vector;
        int nx;

        public SpinResolved_DoubleVector(int nx)
        {
            this.nx = nx;
            spin_vector = new DoubleVector[2];

            spin_vector[0] = new DoubleVector(nx);
            spin_vector[1] = new DoubleVector(nx);
        }

        public double this[int i, Spin spin]
        {
            get { return this.Spin_Vector(spin)[i]; }
            set { this.Spin_Vector(spin)[i] = value; }
        }

        /// <summary>
        /// Cast from DoubleVector to SpinResolved_DoubleMatrix
        /// Assumes equal spin contributions for up and down in mat
        /// </summary>
        public static explicit operator SpinResolved_DoubleVector(DoubleVector vec)
        {
            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(vec.Length);
            result.Spin_Up = vec / 2.0;
            result.Spin_Down = vec / 2.0;

            return result;
        }

        public static SpinResolved_DoubleVector operator *(double scalar, SpinResolved_DoubleVector vec)
        {
            vec.Spin_Up *= scalar; vec.Spin_Down *= scalar;

            return vec;
        }

        public static SpinResolved_DoubleVector operator *(SpinResolved_DoubleVector vec, double scalar)
        {
            return scalar * vec;
        }


        public static SpinResolved_DoubleVector operator +(SpinResolved_DoubleVector vec_1, SpinResolved_DoubleVector vec_2)
        {
            // check if lengths are the same
            if (vec_1.Nx != vec_2.Nx)
                throw new ArrayTypeMismatchException();

            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(vec_1.Nx);
            result.Spin_Up = vec_1.Spin_Up + vec_2.Spin_Up;
            result.Spin_Down = vec_1.Spin_Down + vec_2.Spin_Down;

            return result;
        }

        public static SpinResolved_DoubleVector operator -(SpinResolved_DoubleVector vec_1, SpinResolved_DoubleVector vec_2)
        {
            // check if lengths are the same
            if (vec_1.Nx != vec_2.Nx)
                throw new ArrayTypeMismatchException();

            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(vec_1.Nx);
            result.Spin_Up = vec_1.Spin_Up - vec_2.Spin_Up;
            result.Spin_Down = vec_1.Spin_Down - vec_2.Spin_Down;

            return result;
        }

        /// <summary>
        /// returns the DoubleMatrix for the given spin
        /// </summary>
        public DoubleVector Spin_Vector(Spin spin)
        {
            if (spin == Spin.Up)
                return spin_vector[0];
            else
                return spin_vector[1];
        }

        public DoubleVector Spin_Summed_Vector
        {
            get { return spin_vector[0] + spin_vector[1]; }
        }

        public int Nx
        {
            get { return nx; }
        }

        public DoubleVector Spin_Up
        {
            get { return spin_vector[0]; }
            set { spin_vector[0] = value; }
        }

        public DoubleVector Spin_Down
        {
            get { return spin_vector[1]; }
            set { spin_vector[1] = value; }
        }
    }


    /// <summary>
    /// A class for holding spin resolved DoubleMatrix classes from CenterSpace types
    /// This is just a pair of DoubleMatrices which can be indexed
    /// </summary>
    public class SpinResolved_DoubleMatrix
    {
        int nx, ny;
        DoubleMatrix[] spin_matrix;

        public SpinResolved_DoubleMatrix(int nx, int ny)
        {
            this.nx = nx; this.ny = ny;
            spin_matrix = new DoubleMatrix[2];

            spin_matrix[0] = new DoubleMatrix(nx, ny);
            spin_matrix[1] = new DoubleMatrix(nx, ny);
        }

        public double this[int i, int j, Spin spin]
        {
            get { return this.Spin_Matrix(spin)[i, j]; }
            set { this.Spin_Matrix(spin)[i, j] = value; }
        }

        /// <summary>
        /// Cast from DoubleMatrix to SpinResolved_DoubleMatrix
        /// Assumes equal spin contributions for up and down in mat
        /// </summary>
        public static explicit operator SpinResolved_DoubleMatrix(DoubleMatrix mat)
        {
            SpinResolved_DoubleMatrix result = new SpinResolved_DoubleMatrix(mat.Rows, mat.Cols);
            result.Spin_Up = mat / 2.0;
            result.Spin_Down = mat / 2.0;

            return result;
        }

        public static SpinResolved_DoubleMatrix operator *(double scalar, SpinResolved_DoubleMatrix mat)
        {
            mat.Spin_Up *= scalar; mat.Spin_Down *= scalar;

            return mat;
        }

        public static SpinResolved_DoubleMatrix operator *(SpinResolved_DoubleMatrix mat, double scalar)
        {
            return scalar * mat;
        }

        public static SpinResolved_DoubleMatrix operator +(SpinResolved_DoubleMatrix mat_1, SpinResolved_DoubleMatrix mat_2)
        {
            // check if lengths are the same
            if (mat_1.Nx != mat_2.Nx || mat_1.Ny != mat_2.Ny)
                throw new ArrayTypeMismatchException();

            SpinResolved_DoubleMatrix result = new SpinResolved_DoubleMatrix(mat_1.Nx, mat_1.Ny);
            result.Spin_Up = mat_1.Spin_Up + mat_2.Spin_Up;
            result.Spin_Down = mat_1.Spin_Down + mat_2.Spin_Down;

            return result;
        }

        public static SpinResolved_DoubleMatrix operator -(SpinResolved_DoubleMatrix mat_1, SpinResolved_DoubleMatrix mat_2)
        {
            // check if lengths are the same
            if (mat_1.Nx != mat_2.Nx || mat_1.Ny != mat_2.Ny)
                throw new ArrayTypeMismatchException();

            SpinResolved_DoubleMatrix result = new SpinResolved_DoubleMatrix(mat_1.Nx, mat_1.Ny);
            result.Spin_Up = mat_1.Spin_Up - mat_2.Spin_Up;
            result.Spin_Down = mat_1.Spin_Down - mat_2.Spin_Down;

            return result;
        }

        /// <summary>
        /// returns the DoubleMatrix for the given spin
        /// </summary>
        public DoubleMatrix Spin_Matrix(Spin spin)
        {
            if (spin == Spin.Up)
                return spin_matrix[0];
            else
                return spin_matrix[1];
        }

        public DoubleMatrix Spin_Summed_Matrix
        {
            get { return spin_matrix[0] + spin_matrix[1]; }
        }

        public int Nx
        {
            get { return nx; }
        }

        public int Ny
        {
            get { return ny; }
        }

        public DoubleMatrix Spin_Up
        {
            get { return spin_matrix[0]; }
            set { spin_matrix[0] = value; }
        }

        public DoubleMatrix Spin_Down
        {
            get { return spin_matrix[1]; }
            set { spin_matrix[1] = value; }
        }
    }

    public class SpinResolved_Data
    {
        Band_Data[] spin_data;

        public SpinResolved_Data(Band_Data spin_up, Band_Data spin_down)
        {
            spin_data = new Band_Data[2];
            spin_data[0] = spin_up; spin_data[1] = spin_down;
        }

        public SpinResolved_Data(int nx)
        {
            spin_data = new Band_Data[2];
            spin_data[0] = new Band_Data(nx, 0.0); spin_data[1] = new Band_Data(nx, 0.0);
        }

        public SpinResolved_Data(int nx, int ny)
        {
            spin_data = new Band_Data[2];
            spin_data[0] = new Band_Data(nx, ny, 0.0); spin_data[1] = new Band_Data(nx, ny, 0.0);
        }

        public SpinResolved_Data(int nx, int ny, int nz)
        {
            spin_data = new Band_Data[2];
            spin_data[0] = new Band_Data(nx, ny, nz, 0.0); spin_data[1] = new Band_Data(nx, ny, nz, 0.0);
        }

        /// <summary>
        /// Cast from Band_Data to SpinResolved_Data
        /// Assumes equal spin contributions for up and down
        /// </summary>
        public static explicit operator SpinResolved_Data(Band_Data data)
        {
            SpinResolved_Data result = new SpinResolved_Data(data / 2.0, data / 2.0);

            return result;
        }

        public static SpinResolved_Data operator *(double scalar, SpinResolved_Data data)
        {
            return new SpinResolved_Data(scalar * data.Spin_Up, scalar * data.Spin_Down);
        }

        public static SpinResolved_Data operator /(SpinResolved_Data data, double scalar)
        {
            return (1.0 / scalar) * data;
        }

        public static SpinResolved_Data operator *(SpinResolved_Data data, double scalar)
        {
            return scalar * data;
        }

        public static SpinResolved_Data operator +(SpinResolved_Data data_1, SpinResolved_Data data_2)
        {
            // check if lengths are the same
            if (data_1.Length != data_2.Length || data_1.Spin_Up.Dimension != data_2.Spin_Up.Dimension)
                throw new ArrayTypeMismatchException();

            return new SpinResolved_Data(data_1.Spin_Up + data_2.Spin_Up, data_1.Spin_Down + data_2.Spin_Down);
        }

        public static SpinResolved_Data operator -(SpinResolved_Data data_1, SpinResolved_Data data_2)
        {
            // check if lengths are the same
            if (data_1.Length != data_2.Length || data_1.Spin_Up.Dimension != data_2.Spin_Up.Dimension)
                throw new ArrayTypeMismatchException();

            return new SpinResolved_Data(data_1.Spin_Up - data_2.Spin_Up, data_1.Spin_Down - data_2.Spin_Down);
        }

        public SpinResolved_Data DeepenThisCopy()
        {
            int dim = this.Spin_Down.Dimension;

            if (dim == 1)
            {
                int nz = this.Spin_Down.vec.Length;
                SpinResolved_Data result = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));
                for (int i = 0; i < nz; i++)
                {
                    result.Spin_Down.vec[i] = this.Spin_Down.vec[i];
                    result.Spin_Up.vec[i] = this.Spin_Up.vec[i];
                }

                return result;
            }
            else if (dim == 2)
            {
                int ny = this.Spin_Down.mat.Rows;
                int nz = this.Spin_Down.mat.Cols;
                SpinResolved_Data result = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny, nz)), new Band_Data(new DoubleMatrix(ny, nz)));
                for (int i = 0; i < ny; i++)
                    for (int j = 0; j < nz; j++)
                    {
                        result.Spin_Down.mat[i, j] = this.Spin_Down.mat[i, j];
                        result.Spin_Up.mat[i, j] = this.Spin_Up.mat[i, j];
                    }

                return result;
            }
            else if (dim == 3)
            {
                int nx = this.Spin_Summed_Data.vol[0].Rows;
                int ny = this.Spin_Summed_Data.vol[0].Cols;
                int nz = this.Spin_Summed_Data.vol.Length;

                SpinResolved_Data result = new SpinResolved_Data(new Band_Data(nx, ny, nz, 0.0), new Band_Data(nx, ny, nz, 0.0));

                for (int k = 0; k < nz; k++)
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                        {
                            result.Spin_Up.vol[k][i, j] = this.Spin_Up.vol[k][i, j];
                            result.Spin_Down.vol[k][i, j] = this.Spin_Down.vol[k][i, j];
                        }

                return result;
            }
            else
                throw new NotImplementedException();
        }

        /// <summary>
        /// returns the DoubleMatrix for the given spin
        /// </summary>
        public Band_Data Spin_Data(Spin spin)
        {
            if (spin == Spin.Up)
                return spin_data[0];
            else
                return spin_data[1];
        }

        public Band_Data Spin_Summed_Data
        {
            get { return spin_data[0] + spin_data[1]; }
        }

        public Band_Data Spin_Up
        {
            get { return spin_data[0]; }
            set { spin_data[0] = value; }
        }

        public Band_Data Spin_Down
        {
            get { return spin_data[1]; }
            set { spin_data[1] = value; }
        }

        public int Length
        {
            get { return spin_data[0].Length + spin_data[1].Length; }
        }
    }

}
