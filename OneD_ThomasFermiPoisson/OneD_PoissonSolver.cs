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
        public OneD_PoissonSolver(double dz, int nz, double top_bc, double bottom_bc, double[] layer_depths, bool using_flexPDE, string flexPDE_input, bool freeze_out_dopents, double tol)
            : base (1.0, 1.0, dz, 1, 1, nz, using_flexPDE, flexPDE_input, freeze_out_dopents, tol)
        {
            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // (as phi = E_c / (-1.0 * q_e)
            this.top_bc = top_bc / (-1.0 * Physics_Base.q_e); this.bottom_bc = bottom_bc / (-1.0 * Physics_Base.q_e);

            // generate Laplacian matrix (spin-resolved)
            if (!flexPDE)
            {
                laplacian = Generate_Laplacian();
                lu_fact = new DoubleLUFact(laplacian);
            }
            else
                Create_FlexPDE_Input_File(flexPDE_input, dens_filename, layer_depths);
        }

        public DoubleVector Get_Band_Energy(DoubleVector density)
        {
            if (flexPDE)
                // calculate band energy using a potential found by calling FlexPDE
                return Get_BandEnergy_From_FlexPDE(new Band_Data(density), dens_filename).vec;
            else
                // calculate band energies using a potential calculated on a regular grid (not ideal, or scalable)
                return Get_BandEnergy_On_Regular_Grid(density);
        }

        protected override void Save_Density(Band_Data density, string filename)
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

        protected override Band_Data Parse_Potential(string[] data, int first_line)
        {
            // and check that there is the right number of data points back
            if (data.Length - first_line != nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleVector(nz));
            for (int i = 0; i < nz; i++)
                result.vec[i] = double.Parse(data[first_line + i]);

            return result;
        }

        DoubleVector Get_BandEnergy_On_Regular_Grid(DoubleVector density)
        {
            // set the top and bottom boundary conditions
            double factor = -1.0 * Physics_Base.epsilon / (dz * dz);
            density[0] = top_bc * factor; density[density.Length - 1] = bottom_bc * factor;

            // solve Poisson's equation
            DoubleVector potential = lu_fact.Solve(density);
            
            // return band energy
            return -1.0 * Physics_Base.q_e * potential;
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
        protected override void Create_FlexPDE_Input_File(string flexPDE_input, string dens_filename, double[] layer_depths)
        {
            string well_dens_filename = dens_filename;
            string dopent_dens_filename = dens_filename;

            double dopent_depth = layer_depths[1];
            double dopent_width = layer_depths[2] - layer_depths[1];
            double well_depth = layer_depths[layer_depths.Length - 2];           // put the top of the heterostructure boundary as the penultimate layer boundary

            // check if the dopents should be frozen out
            if (freeze_out_dopents)
                // if true, change the input density filename for the dopents
                dopent_dens_filename = "dopents_frozen_" + dens_filename;
            
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
            // this is where the density variable
            sw.WriteLine("\trho\t! density");
            // number of lattice sites that the density needs to be output to
            sw.WriteLine("\tnx = " + nz.ToString());
            // size of the sample
            sw.WriteLine("\tlx = " + (nz * dz).ToString());
            // device dimensions
            sw.WriteLine("\tdopent_depth = " + dopent_depth.ToString());
            sw.WriteLine("\tdopent_width = " + dopent_width.ToString());
            sw.WriteLine("\twell_depth = " + well_depth.ToString());
            // the top boundary condition on the surface of the sample
            sw.WriteLine("\t! Boundary conditions");
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine("\tbottom_V = " + bottom_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivities");
            sw.WriteLine("\teps_0 = " + Physics_Base.epsilon_0.ToString());
            // relative permitivity of GaAs
            sw.WriteLine("\teps_r = " + Physics_Base.epsilon_r.ToString());
            sw.WriteLine("\teps");
            sw.WriteLine();
            // boundary layer definitions (note that this is quite a specific layer structure)
            sw.WriteLine("\t! boundary layer definitions");
	        sw.WriteLine("\tdopent_top = " + dopent_depth.ToString());
	        sw.WriteLine("\tdopent_bottom = " + (dopent_depth + dopent_width).ToString());
            sw.WriteLine("\twell_top = " + well_depth.ToString());           // put the top of the heterostructure boundary as the penultimate layer boundary
            sw.WriteLine("\twell_bottom = lx");
            sw.WriteLine(); 
            sw.WriteLine("EQUATIONS");
            // Poisson's equation (not too happy about this... shouldn't it be -1.0 * rho?!)
            sw.WriteLine("\tu: div(eps * grad(u)) = rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");
            sw.WriteLine("\tREGION 1	! capping layer");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0 * eps_r");
            sw.WriteLine("\t\tSTART(0)");
            sw.WriteLine("\t\tPOINT VALUE(u) = top_V");
            sw.WriteLine("\t\tLINE TO (dopent_top)");
            sw.WriteLine("\tREGION 2	! dopent layer");
            sw.WriteLine("\t\trho = TABLE(\'" + dopent_dens_filename + "\', x)");
            sw.WriteLine("\t\teps = eps_0 * eps_r");
            sw.WriteLine("\t\tSTART(dopent_top)");
            sw.WriteLine("\t\tLINE TO (dopent_bottom)");
            sw.WriteLine("\tREGION 3	! separator layer");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0 * eps_r");
            sw.WriteLine("\t\tSTART(dopent_bottom)");
            sw.WriteLine("\t\tLINE TO (well_top)");
            sw.WriteLine("\tREGION 3	! well layer");
            sw.WriteLine("\t\trho = TABLE(\'" + well_dens_filename + "\', x)");
            sw.WriteLine("\t\teps = eps_0 * eps_r");
            sw.WriteLine("\t\tSTART(well_top)");
            sw.WriteLine("\t\tLINE TO (well_bottom)");
            // a possible bottomr boundary condition which hasn't worked yet...
            // this form of the boundary condition is equivalent to " d(eps * u)/dn + (bottom_V - u) = 0.0 "
            // which gives a combined Neumann/Dirichlet 
            sw.WriteLine("\t\tPOINT VALUE(u) = bottom_V");
            //sw.WriteLine("\t\tPOINT NATURAL(u) = (bottom_V - u)");
            sw.WriteLine();
            sw.WriteLine("PLOTS");
            sw.WriteLine("\tELEVATION(rho) FROM (0) TO (lx)");
	        sw.WriteLine("\tELEVATION(u) FROM (0) TO (lx) EXPORT(nx) FORMAT \'#1\' FILE=\'pot.dat\'");
            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }
    }
}
