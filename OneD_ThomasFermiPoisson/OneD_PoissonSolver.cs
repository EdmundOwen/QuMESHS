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
        double top_bc, bottom_bc;
        // this is where the density will be saved out to
        string dens_filename = "dens.dat";

        // parameters for regular grid solve
        DoubleMatrix laplacian;
        DoubleLUFact lu_fact;

        // 
        public OneD_PoissonSolver(double dz, int nz, double top_bc, double bottom_bc, bool using_flexPDE, string flexPDE_input, double tol)
            : base (1.0, 1.0, dz, 1, 1, nz, using_flexPDE, flexPDE_input, tol)
        {
            this.top_bc = top_bc; this.bottom_bc = bottom_bc;

            // generate Laplacian matrix (spin-resolved)
            if (!flexPDE)
            {
                laplacian = Generate_Laplacian();
                lu_fact = new DoubleLUFact(laplacian);
            }
            else
                Create_FlexPDE_Input_File(flexPDE_input, dens_filename);
        }

        public DoubleVector Get_Potential(DoubleVector density)
        {
            if (flexPDE)
                // calculate potential by calling FlexPDE
                return Get_Potential_From_FlexPDE(new Potential_Data(density), dens_filename).vec;
            else
                // calculate potential on a regular grid (not ideal, or scalable)
                return Get_Potential_On_Regular_Grid(density);
        }

        protected override void Save_Density(Potential_Data density, string filename)
        {
            // check that the dimension of the density is correct
            if (density.Dimension != 1)
                throw new RankException();

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
                sw.Write(((float)density.vec[i]).ToString() + '\t');

            sw.Close();
        }

        protected override Potential_Data Parse_Potential(string[] data, int first_line)
        {
            // and check that there is the right number of data points back
            if (data.Length - first_line != nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Potential_Data result = new Potential_Data(new DoubleVector(nz));
            for (int i = 0; i < nz; i++)
                result.vec[i] = double.Parse(data[first_line + i]);

            return result;
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
            double factor = -1.0 * Physics_Base.epsilon / (dz * dz);

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

        /// <summary>
        /// creates an input file for flexPDE to solve a 1D poisson equation
        /// </summary>
        protected override void Create_FlexPDE_Input_File(string flexPDE_input, string dens_filename)
        {
            // check if an input file already exists and delete it
            if (File.Exists(flexPDE_input))
                File.Delete(flexPDE_input);

            // open up a new streamwriter to create the input file
            StreamWriter sw = new StreamWriter(flexPDE_input);

            // write the file
            sw.WriteLine("TITLE \'Band Structure\'");
            sw.WriteLine("COORDINATES cartesian1");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            sw.WriteLine("\tERRLIM=1e-5");
            sw.WriteLine("DEFINITIONS");
            // this is where the input data for the density is defined
            sw.WriteLine("\trho = TABLE(\'" + dens_filename + "\', x)");
            // number of lattice sites that the density needs to be output to
            sw.WriteLine("\tnx = " + nz.ToString());
            // size of the sample
            sw.WriteLine("\tlx = " + (nz * dz).ToString());
            // the top boundary condition on the surface of the sample
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine("\tbottom_V = " + bottom_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\teps_0 = 1.41859713");
            // relative permitivity of GaAs
            sw.WriteLine("\teps_r = 13");
            sw.WriteLine("\teps");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: div(eps * grad(u)) = rho");
            sw.WriteLine();
            sw.WriteLine("BOUNDARIES");
            sw.WriteLine("\tREGION 1");
            sw.WriteLine("\t\teps = eps_0 * eps_r");
            sw.WriteLine("\t\tSTART(0)");
            sw.WriteLine("\t\tPOINT VALUE(u) = top_V");
            sw.WriteLine("\t\tLINE TO (lx)");
            // this form of the boundary condition is equivalent to " 1.0 * d(eps * u)/dn + (1.0 / bottom_v) * u = 0.0 "
            // which gives a combined Neumann/Dirichlet 
            sw.WriteLine("\t\tPOINT NATURAL(u) = (1.0 / bottom_V) * u");
            sw.WriteLine();
            sw.WriteLine("PLOTS");
            sw.WriteLine("\tELEVATION(rho) FROM (0) TO (lx)");
	        sw.WriteLine("\tELEVATION(u) FROM (0) TO (lx) export(nx) format \'#1\' file=\'pot.dat\'");
            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }
    }
}
