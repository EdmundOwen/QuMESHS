using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;

namespace OneD_ThomasFermiPoisson
{
    enum RBF_Type
    {
        Gaussian,
    }

    /// <summary>
    /// A class for calculating a one dimensional radial basis function fit for interpolation purposes
    /// </summary>
    public class OneD_RBF_Fit
    {
        double[] positions;
        double width;

        DoubleVector weights;
        DoubleMatrix svd_mat;
        DoubleSVDLeastSq svd_solv;
        RBF_Type rbf_type = RBF_Type.Gaussian;

        public OneD_RBF_Fit(double[] positions, double width)
        {
            this.positions = positions; this.width = width;
        }

        public OneD_RBF_Fit(double[] donor_dens, double[] positions, double width)
        {
            this.positions = positions; this.width = width;

            // generate the svd matrix
            svd_mat = Generate_SVD_Matrix(positions);
            // and solve it to get the weights (this is simply solving A x = y where x are the weights and y are the data values)
            svd_solv = new DoubleSVDLeastSq(svd_mat);
            weights = svd_solv.Solve(new DoubleVector(donor_dens));
        }

        /// <summary>
        /// generate a matrix of radial basis functions
        /// </summary>
        DoubleMatrix Generate_SVD_Matrix(double[] positions)
        {
            int order = positions.Length;
            DoubleMatrix result = new DoubleMatrix(order, order);

            for (int i = 0; i < order; i++)
                for (int j = 0; j < order; j++)
                    result[i, j] = Get_RBF(positions[i], positions[j]);

            return result;
        }

        /// <summary>
        /// get radial basis function value at a given distance
        /// </summary>
        double Get_RBF(double pos_diff)
        {
            switch (rbf_type)
            {
                case RBF_Type.Gaussian:
                    return Math.Exp(-1.0 * pos_diff * pos_diff / width);

                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// gets the radial basis function value for a given pair of positions
        /// </summary>
        double Get_RBF(double pos_1, double pos_2)
        {
            return Get_RBF(Math.Abs(pos_1 - pos_2));
        }

        /// <summary>
        /// returns the RBF weights for the given fitting values
        /// </summary>
        public DoubleVector Get_RBF_Weights(double[] fit_vals)
        {
            svd_mat = Generate_SVD_Matrix(positions);
            svd_solv = new DoubleSVDLeastSq(svd_mat);

            weights = svd_solv.Solve(new DoubleVector(fit_vals));
            return weights;
        }

        /// <summary>
        /// gets a string with rbf equation
        /// </summary>
        /// <returns></returns>
        public string Get_RBF_Equation(string output_variable, string input_variable)
        {
            if (weights == null)
                throw new Exception("Error - you have not calculated the weights yet!");
            
            string result = output_variable + " = ";

            for (int i = 0; i < weights.Length - 1; i++)
            {
                result += weights[i].ToString() + " * " + Get_RBF_String(input_variable, i);
                if (weights[i + 1] > 0.0)
                    result += " + ";
            }

            // and add the final value
            result += weights[weights.Length].ToString() + " * " + Get_RBF_String(input_variable, weights.Length);

            return result;
        }

        /// <summary>
        /// gets a string with the given RBF number
        /// </summary>
        string Get_RBF_String(string input_variable, int rbf_num)
        {
            switch (rbf_type)
            {
                case RBF_Type.Gaussian:
                    return "exp(-1.0 * (" + positions[rbf_num].ToString() + " - " + input_variable + ")^2 / " + width.ToString() + ")";

                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// prints the weights to the specified file
        /// </summary>
        public void Print_Weights(string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            for (int i = 0; i < weights.Length; i++)
                sw.WriteLine(weights[i].ToString());

            sw.Close();
        }

        /// <summary>
        /// width of the radial basis function
        /// </summary>
        public double Width
        {
            get { return width; }
            set
            {
                if (value <= 0.0) throw new Exception("Error - widths of RBFs must be positive!");
                else width = value;
            }
        }
    }
}
