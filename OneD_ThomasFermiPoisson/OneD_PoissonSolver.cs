using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using Solver_Bases;
using System.Threading;

namespace OneD_ThomasFermiPoisson
{
    class OneD_PoissonSolver : Potential_Base
    {
        // parameters for regular grid solve
        int nz;
        double dz;

        double top_bc, bottom_bc;

        bool flexPDE;
        string flexpde_inputfile;

        DoubleMatrix laplacian;
        DoubleLUFact lu_fact;

        // 
        public OneD_PoissonSolver(double dz, int nz, double top_bc, double bottom_bc, bool using_flexPDE, string flexPDE_input)
        {
            this.nz = nz; this.dz = dz;

            this.top_bc = top_bc; this.bottom_bc = bottom_bc;

            // check whether using flexPDE
            flexPDE = using_flexPDE;
            if (using_flexPDE)
                this.flexpde_inputfile = flexPDE_input;

            // generate Laplacian matrix (spin-resolved)
            if (!flexPDE)
            {
                laplacian = Generate_Laplacian();
                lu_fact = new DoubleLUFact(laplacian);
            }
        }

        public DoubleVector Get_Potential(DoubleVector density)
        {
            if (flexPDE)
                // calculate potential by calling FlexPDE
                return Get_Potential_From_FlexPDE(density);
            else
                // calculate potential on a regular grid (not ideal, or scalable)
                return Get_Potential_On_Regular_Grid(density);
        }

        /// <summary>
        /// gets the potential using flexPDE
        /// </summary>
        DoubleVector Get_Potential_From_FlexPDE(DoubleVector density)
        {
            // save density to file in a FlexPDE "TABLE" format
            Save_Density(density, "density_1d.dat");
            // remove pot.dat if it still exists (to make sure that a new data file is made by flexPDE)
            try { File.Delete("pot.dat"); } catch (Exception) { }

            // run the flexPDE program as a process (quietly)
            Process.Start("C:\\FlexPDE6\\FlexPDE6.exe", "-Q " + flexpde_inputfile);
            while (!File.Exists("pot.dat"))
                Thread.Sleep(10);
            Thread.Sleep(10);
            
            string[] lines = File.ReadAllLines("pot.dat");
            // work out where the data starts (this is flexPDE specific)
            int first_line = 0;
            for (int i = 0; i < lines.Length; i++)
                if (lines[i].StartsWith("}"))
                {
                    first_line = i + 1;
                    break;
                }

            // and check that there is the right number of data points back
            if (lines.Length - first_line != nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            DoubleVector result = new DoubleVector(nz);
            for (int i = 0; i < nz; i++)
                result[i] = double.Parse(lines[first_line + i]);

            return result;
        }

        void Save_Density(DoubleVector density, string filename)
        {
            // open stream
            StreamWriter sw = new StreamWriter(filename);

            // save out positions
            sw.WriteLine("x " + nz.ToString());
            for (int i = 0; i < nz;i++)
                sw.Write(((float)(i * dz)).ToString() + '\t');

            // save out densities
            sw.WriteLine();
            sw.WriteLine("data");
            for (int i = 0; i < nz; i++)
                sw.Write(((float)density[i]).ToString() + '\t');

            sw.Close();
        }

        DoubleVector Get_Potential_On_Regular_Grid(DoubleVector spin_resolved_density)
        {
            // sum the spin contributions of the density
            DoubleVector density = new DoubleVector(nz);
            for (int i = 0; i < nz; i++)
                density[i] = spin_resolved_density[i] + spin_resolved_density[i + nz];

            // set boundary conditions
            density[0] = top_bc; density[nz - 1] = bottom_bc;

            DoubleVector potential = lu_fact.Solve(density);
            
            return potential;
        }

        /// <summary>
        /// Generates a spin-resolved laplacian matrix in one-dimension on a regular grid with Dirichlet BCs
        /// </summary>
        DoubleMatrix Generate_Laplacian()
        {
            // the factor which multiplies the Laplace equation
            double factor = -1.0 * epsilon / (dz * dz);

            DoubleMatrix result = new DoubleMatrix(nz, nz);
            for (int i = 0; i < nz - 1; i++)
            {
                // on-diagonal term
                result[i, i] = 2.0 * factor;
                // off-diagonal
                result[i + 1, i] = -1.0 * factor;
                result[i, i + 1] = -1.0 * factor;
            }

            // and fix boundary conditions
            result[0, 0] = 1.0 * factor;
            result[0, 1] = 0.0;
            result[nz - 1, nz - 1] = 1.0 * factor;
            result[nz - 1, nz - 2] = 0.0;

            return result;
        }
    }
}
