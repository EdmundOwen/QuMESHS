using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_PoissonSolver : Potential_Base
    {
        public TwoD_PoissonSolver(double dy, double dz, int ny, int nz, bool using_flexPDE, string flexPDE_input)
            : base(1.0, dy, dz, 1, ny, nz, using_flexPDE, flexPDE_input)
        {
        }

        public DoubleMatrix Get_Potential(DoubleMatrix density)
        {
            if (flexPDE)
                // calculate potential by calling FlexPDE
                return Get_Potential_From_FlexPDE(new Potential_Data(density)).mat;
            else
                // calculate potential on a regular grid (not ideal, or scalable)
                throw new NotImplementedException();// return Get_Potential_On_Regular_Grid(density);
        }

        protected override void Save_Density(Potential_Data density, string filename)
        {
            // check that the dimension of the density is correct
            if (density.Dimension != 2)
                throw new RankException();

            // open stream
            StreamWriter sw = new StreamWriter(filename);

            // save out positions
            sw.WriteLine("x " + nz.ToString());
            for (int i = 0; i < ny; i++)
                //for (int j = 0; j < nz; j++)
                    sw.Write(((float)(i * dy)).ToString() + '\t');

            sw.WriteLine();
            sw.WriteLine();
            sw.WriteLine("y " + nz.ToString());
            //for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    sw.Write(((float)(j * dz)).ToString() + '\t');

            // save out densities
            sw.WriteLine();
            sw.WriteLine("data");
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    sw.Write(((float)density.mat[j, i]).ToString() + '\t');

            sw.Close();
        }

        protected override Potential_Data Parse_Potential(string[] data, int first_line)
        {
            // and check that there is the right number of data points back
            if (data.Length - first_line != ny * nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Potential_Data result = new Potential_Data(new DoubleMatrix(ny, nz));
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    result.mat[i, j] = double.Parse(data[first_line + j * ny + i]);

            return result;
        }

        /*public DoubleVector Get_Well_Potential(int nx)
        {
            Console.WriteLine("THIS IS A PROTOTYPE! NO POISSON SOLVER IMPLEMENTED!!!");

            DoubleVector result = new DoubleVector(nx);
            for (int i = 0; i < nx; i++)
                result[i] = (i - (double)(nx - 1) / 2) * (i - (double)(nx - 1) / 2);

            return result;
        }*/
    }
}
