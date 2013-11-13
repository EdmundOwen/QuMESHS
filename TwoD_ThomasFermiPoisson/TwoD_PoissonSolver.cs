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
        string dens_filename = "dens_2D.dat";
        double bottom_bc;

        public TwoD_PoissonSolver(double dy, double dz, int ny, int nz, double bottom_bc, double[] layer_depths, bool using_flexPDE, string flexPDE_input, string flexPDE_location, bool freeze_out_dopents, double tol)
            : base(1.0, dy, dz, 1, ny, nz, using_flexPDE, flexPDE_input, flexPDE_location, freeze_out_dopents, tol)
        {
            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // (as phi = E_c / (-1.0 * q_e))
            this.bottom_bc = bottom_bc / (-1.0 * Physics_Base.q_e);

            Create_FlexPDE_Input_File(flexPDE_input, dens_filename, layer_depths);
        }

        public DoubleMatrix Get_Band_Energy(DoubleMatrix density)
        {
            if (flexPDE)
                // calculate potential by calling FlexPDE
                return Get_BandEnergy_From_FlexPDE(new Band_Data(density), dens_filename).mat;
            else
                // calculate potential on a regular grid (not ideal, or scalable)
                throw new NotImplementedException();
        }

        public override void Save_Density(Band_Data density, string filename)
        {
            // check that the dimension of the density is correct
            if (density.Dimension != 2)
                throw new RankException();

            // open stream
            StreamWriter sw = new StreamWriter(filename);

            // save out positions
            sw.WriteLine("x " + ny.ToString());
            for (int i = 0; i < ny; i++)
                //for (int j = 0; j < nz; j++)
                sw.Write(((float)(i - ny / 2) * dy).ToString() + '\t');

            sw.WriteLine();
            sw.WriteLine();
            sw.WriteLine("y " + nz.ToString());
            for (int j = nz - 1; j >= 0; j--)
            //for (int j = 0; j < nz; j++)
                sw.Write(((float)(-j * dz)).ToString() + '\t');

            // save out densities
            sw.WriteLine();
            sw.WriteLine("data");
            //for (int i = 0; i < nz; i++)
            for (int i = nz - 1; i >= 0; i--)
            {
                for (int j = 0; j < ny; j++)
                    if (Math.Abs(density.mat[j, i]) < 1e-20)
                        sw.Write("0\t");
                    else
                        // note that the ordering is y first, then z -- this is FlexPDE specific
                        sw.Write(((float)density.mat[j, i]).ToString() + '\t');
                sw.WriteLine();
            }

            sw.Close();
        }

        protected override Band_Data Parse_Potential(string[] data, int first_line)
        {
            // and check that there is the right number of data points back
            if (data.Length - first_line != ny * nz)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleMatrix(ny, nz));
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    result.mat[i, j] = double.Parse(data[first_line + j * ny + i]);

            return result;
        }

        
        /// <summary>
        /// creates an input file for flexPDE to solve a 2D poisson equation (not yet implemented)
        /// </summary>
        protected override void Create_FlexPDE_Input_File(string flexPDE_input, string dens_filename, double[] layer_depths)
        {
            string well_dens_filename = dens_filename;
            string dopent_dens_filename = dens_filename;
            
            double pmma_depth = 200;
            double gate_depth = 10;
            double top_depth = 10;
  
            double split_width = 700;

            // as with OneD_PoissonSolver, this assumes that the layer structure is
            // cap --> dopents --> spacer --> heterostructure
            // well_width is ~50 nm where there is substantial density
            double dopent_depth = layer_depths[1];
            double dopent_width = layer_depths[2] - layer_depths[1];
            double well_depth = layer_depths[layer_depths.Length - 2];
            double well_width = 50;

            double split_V = -0.0;
            double top_V = -0.0;

            // check if the dopents should be frozen out
            if (freeze_out_dopents)
                // if true, change the input density filename for the dopents
                dopent_dens_filename = "dopents_frozen_" + dens_filename;

            // check if an input file already exists and delete it
            if (File.Exists(flexPDE_input))
                File.Delete(flexPDE_input);

            // open up a new streamwriter to create the input file
            StreamWriter sw = new StreamWriter(flexPDE_input);
            
            sw.WriteLine("TITLE \'Split Gate\'");
            sw.WriteLine("COORDINATES cartesian2");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            //sw.WriteLine("\tERRLIM=1e-6");
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\trho");
            sw.WriteLine();
            // lattice dimensions for outputting the potential for the density solver
            sw.WriteLine("\tny = " + ny.ToString());
            sw.WriteLine("\tnz = " + nz.ToString());
            // simulation dimension
            sw.WriteLine("\tly = " + (dy * ny).ToString());
            sw.WriteLine("\tlz = " + (dz * nz).ToString());
            sw.WriteLine();
            // device dimensions
            sw.WriteLine("\tpmma_depth = " + pmma_depth.ToString());
            sw.WriteLine("\tgate_depth = " + gate_depth.ToString());
            sw.WriteLine("\tsplit_width = " + split_width.ToString());
            sw.WriteLine("\ttop_depth = " + top_depth.ToString());
            sw.WriteLine("\twell_depth = " + well_depth.ToString());
            sw.WriteLine("\twell_width = " + well_width.ToString());
            sw.WriteLine("\tdopent_depth = " + dopent_depth.ToString());
            sw.WriteLine("\tdopent_width = " + dopent_width.ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tsplit_V = " + split_V.ToString());
            sw.WriteLine("\ttop_V = " + top_V.ToString());
            sw.WriteLine();
            // bottom boundary condition
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivities");
            sw.WriteLine("\teps_0 = " + Physics_Base.epsilon_0.ToString());
            // relative permitivity of GaAs
            sw.WriteLine("\teps_r = " + Physics_Base.epsilon_r.ToString());
            sw.WriteLine("\teps_pmma = " + Physics_Base.epsilon_pmma.ToString());
            sw.WriteLine("\teps");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation (not too happy about this... shouldn't it be -1.0 * rho?!)
            sw.WriteLine("\tu: div(eps * grad(u)) = rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");
            sw.WriteLine("\tREGION 1 ! Total simulation domain");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(-ly / 2, -lz)");
            // bottom of device boundary condition
            sw.WriteLine("\t\tVALUE(u) = bottom_bc");
            sw.WriteLine("\t\tLINE TO (ly / 2, -lz)");
            // and reset to default boundary conditions for all other edges
            sw.WriteLine("\t\tNATURAL(u) = 0");
            sw.WriteLine("\t\tLINE TO (ly / 2, pmma_depth + top_depth)");
            sw.WriteLine("\t\tLINE TO (-ly / 2, pmma_depth + top_depth)");
            sw.WriteLine("\t\tLINE TO CLOSE");
            sw.WriteLine("\tREGION 2 ! Material");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0 * eps_r");
            sw.WriteLine("\t\tSTART (-ly / 2, -lz)");
	        sw.WriteLine("\t\tLINE TO (-ly / 2, 0) TO (ly / 2, 0) TO (ly / 2, -lz) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION 3 ! PMMA Layer");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0 * eps_pmma");
	        sw.WriteLine("\t\tSTART (-ly / 2, 0)");
            sw.WriteLine("\t\tLINE TO (-ly / 2, pmma_depth) TO (ly / 2, pmma_depth) TO (ly / 2, 0) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION 4 ! Left split gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0");
	        sw.WriteLine("\t\tSTART(-ly / 2, 0)");
	        // left split gate voltate
            sw.WriteLine("\t\tVALUE(u) = split_V");
	        sw.WriteLine("\t\tLINE TO (-ly / 2, gate_depth) TO (-split_width / 2, gate_depth) TO (-split_width / 2, 0) TO CLOSE");
	        sw.WriteLine();
            sw.WriteLine("\tREGION 5 ! Right split gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(ly / 2, 0)");
            // right split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
	        sw.WriteLine("\t\tLINE TO (ly / 2, gate_depth) TO (split_width / 2, gate_depth) TO (split_width / 2, 0) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION 6 ! Top gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\teps = eps_0");
	        sw.WriteLine("\t\tSTART(ly / 2, pmma_depth)");
            // top gate voltage
	        sw.WriteLine("\t\tVALUE(u) = top_V");
	        sw.WriteLine("\t\tLINE TO (ly / 2, pmma_depth + top_depth) TO (-ly / 2, pmma_depth + top_depth) TO (-ly / 2, pmma_depth) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION 7 ! Quantum well");
            sw.WriteLine("\t\trho = TABLE(\'" + well_dens_filename + "\', x, y)");
	        sw.WriteLine("\t\teps = eps_0 * eps_r");
	        sw.WriteLine("\t\tSTART(ly / 2, -well_depth)");
	        sw.WriteLine("\t\tLINE TO (ly / 2, -lz) TO (-ly / 2, -lz) TO (-ly / 2, -well_depth) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION 8 ! Dopent layer");
            sw.WriteLine("\t\trho = TABLE(\'" + dopent_dens_filename + "\', x, y)");
	        sw.WriteLine("\t\teps = eps_0 * eps_r");
	        sw.WriteLine("\t\tSTART(ly / 2, -dopent_depth)");
	        sw.WriteLine("\t\tLINE TO (ly / 2, -dopent_depth - dopent_width) TO (-ly / 2, -dopent_depth - dopent_width) TO (-ly / 2, -dopent_depth) TO CLOSE");
            sw.WriteLine();
            //{REGION 8 {Substrate}
	        //eps = eps_0 * eps_r
	        //rho = dopent
	        //START(-lx/2, -ly)
	        //LINE TO (-lx/2, -ly-2000) 
            //TO (lx/2, -ly-2000) 
            //TO (lx/2, -ly) TO CLOSE}
            sw.WriteLine();
            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tCONTOUR(u)");
            sw.WriteLine("\tCONTOUR(u) ZOOM (-ly / 2, -lz, ly, lz) EXPORT(ny, nz) FORMAT \'#1\' FILE = \'pot.dat\'");
	        sw.WriteLine("\tELEVATION(u) FROM (-ly / 2, -well_depth - well_width / 2) TO (ly / 2, -well_depth - well_width / 2)");
            sw.WriteLine("\tELEVATION(u) FROM (0, 0) TO (0, -lz)");
            sw.WriteLine("\tELEVATION(rho) FROM (0, 0) TO (0, -lz)");
            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }
    }
}
