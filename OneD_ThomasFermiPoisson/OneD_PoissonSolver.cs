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
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace OneD_ThomasFermiPoisson
{
    public class OneD_PoissonSolver : Potential_Base
    {
        double top_bc, bottom_bc;
        double top_eps, bottom_eps;

        // parameters for regular grid solve
        DoubleMatrix laplacian;
        DoubleLUFact lu_fact;

        /*
        public OneD_PoissonSolver(double dz, int nz, double top_bc, double bottom_bc, double[] layer_depths, bool using_flexPDE, string flexPDE_input, string flexPDE_location, bool freeze_out_dopents, double tol)
            : base (1.0, 1.0, dz, 1, 1, nz, using_flexPDE, flexPDE_input, flexPDE_location, freeze_out_dopents, tol)
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

            this.dens_filename = "dens_1D.dat";

            Create_FlexPDE_Input_File(flexPDE_input, dens_filename, layer_depths);
        }
        */

        public OneD_PoissonSolver(double dz, int nz, ILayer[] layers, bool using_flexPDE, string flexPDE_input, string flexPDE_location, double tol)
            : base(1.0, 1.0, dz, 1, 1, nz, layers, using_flexPDE, flexPDE_input, flexPDE_location, tol)
        {
            // generate Laplacian matrix (spin-resolved)
            if (!flexPDE)
            {
                laplacian = Generate_Laplacian(layers);
                lu_fact = new DoubleLUFact(laplacian);
            }

            this.dens_filename = "dens_1D.dat";
        }

        /*
        public DoubleVector Get_Band_Energy(DoubleVector density)
        {
            if (flexPDE)
                // calculate band energy using a potential found by calling FlexPDE
                return Get_BandEnergy_From_FlexPDE(new Band_Data(density), dens_filename).vec;
            else
                // calculate band energies using a potential calculated on a regular grid (not ideal, or scalable)
                return Get_BandEnergy_On_Regular_Grid(density);
        }
        */

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

        protected override Band_Data Get_BandEnergy_On_Regular_Grid(Band_Data charge_density)
        {
            // set the top and bottom boundary conditions where [0] is the bottom of the device
            charge_density.vec[0] = bottom_bc * -1.0 * bottom_eps / (dz * dz);
            charge_density.vec[charge_density.Length - 1] = top_bc * -1.0 * top_eps / (dz * dz);

            // solve Poisson's equation
            Band_Data potential = new Band_Data(lu_fact.Solve(charge_density.vec));

            // return band energy
            return -1.0 * Physics_Base.q_e * potential;
        }

        /// <summary>
        /// Generates a laplacian matrix in one-dimension on a regular grid with Dirichlet BCs and assuming constant permitivity
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
        /// Generates a laplacian matrix in one-dimension on a regular grid with Dirichlet BCs wot varying permitivity
        /// </summary>
        DoubleMatrix Generate_Laplacian(ILayer[] layers)
        {
            DoubleMatrix result = new DoubleMatrix(nz, nz);
            double factor_plus, factor_minus;

            // cycle through the structure and fill in the Laplacian with the correct permittivities
            double zmin = Geom_Tool.Get_Zmin(layers);
            for (int i = 1; i < nz - 1; i++)
            {
                double eps_plus = Geom_Tool.GetLayer(layers, i * dz + 0.5 * dz + zmin).Permitivity;
                double eps_minus = Geom_Tool.GetLayer(layers, i * dz - 0.5 * dz + zmin).Permitivity;

                // the factor which multiplies the Laplace equation
                factor_plus = -1.0 * eps_plus / (dz * dz);
                factor_minus = -1.0 * eps_minus / (dz * dz);

                // on-diagonal term
                result[i, i] = factor_minus + factor_plus;
                // off-diagonal
                result[i, i - 1] = -1.0 * factor_minus;
                result[i, i + 1] = -1.0 * factor_plus;
            }

            // and fix boundary conditions
            double factor = -1.0 * Geom_Tool.GetLayer(layers, zmin).Permitivity / (dz * dz);
            result[0, 0] = 1.0 * factor;
            result[0, 1] = 0.0;
            factor = -1.0 * Geom_Tool.GetLayer(layers, (nz - 1) * dz + zmin).Permitivity / (dz * dz);
            result[nz - 1, nz - 1] = 1.0 * factor;
            result[nz - 1, nz - 2] = 0.0;

            return result;
        }

        public double Get_Surface_Charge(Band_Data band_offset, ILayer[] layers)
        {
            // calculate the electric field just below the surface
            int surface = (int)(-1.0 * Math.Floor(Geom_Tool.Get_Zmin(layers) / dz));
            double eps = layers[Geom_Tool.Find_Layer_Below_Surface(layers)].Permitivity;
            // by Gauss' theorem, rho = - epsilon_0 * epsilon_r * dV/dz
            double surface_charge = eps * (band_offset[surface] - band_offset[surface - 1]) / dz;
            // divide by q_e to convert the band energy into a potential
            surface_charge /= Physics_Base.q_e;
            // and divide by dz to give a density
            surface_charge /= dz;

            return surface_charge;
        }

        /*
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
        */

        protected override void Create_FlexPDE_Input_File(string flexPDE_input, ILayer[] layers)
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
            // this is where the density variable
            sw.WriteLine("\trho\t! density");
            // number of lattice sites that the density needs to be output to
            sw.WriteLine("\tnx = " + nz.ToString());
            // size of the sample
            sw.WriteLine("\tlx = " + (nz * dz).ToString());
            // the top boundary condition on the surface of the sample
            sw.WriteLine("\t! Boundary conditions");
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine("\tbottom_V = " + bottom_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivities");
            sw.WriteLine("\teps_0 = " + Physics_Base.epsilon_0.ToString());
            // relative permitivity of materials
            sw.WriteLine("\teps_r = " + Physics_Base.epsilon_r.ToString());
            sw.WriteLine("\teps_PMMA = " + Physics_Base.epsilon_pmma.ToString());
            sw.WriteLine("\teps");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation (not too happy about this... shouldn't it be -1.0 * rho?!)
            sw.WriteLine("\tu: div(eps * grad(u)) = rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");

            // cycle through layers
            for (int i = 1; i < layers.Length; i++)
            {
                sw.WriteLine("\tREGION " + layers[i].Layer_No.ToString());
                sw.WriteLine("\t\trho = TABLE(\'" + dens_filename + "\', x)");
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(layers[i].Material));
                sw.WriteLine("\t\tSTART(" + layers[i].Zmin.ToString() + ")");
                if (i == 1)
                    sw.WriteLine("\t\tPOINT VALUE(u) = top_V");
                sw.WriteLine("\t\tLINE TO (" + layers[i].Zmax.ToString() + ")");
                if (i == layers.Length - 1)
                    sw.WriteLine("\t\tPOINT VALUE(u) = bottom_V");
                sw.WriteLine();
            }

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tELEVATION(rho) FROM (" + zmin.ToString() + ") TO (" + (zmin + nz * dz).ToString() + ")");
            sw.WriteLine("\tELEVATION(u) FROM (" + zmin.ToString() + ") TO (" + (zmin + nz * dz).ToString() + ") EXPORT(nx) FORMAT \'#1\' FILE=\'pot.dat\'");
            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }

        public void Set_Boundary_Conditions(ILayer[] layers, double top_bc, double bottom_bc, double top_pos, double bottom_pos)
        {
            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // (as phi = E_c / (-1.0 * q_e)
            this.top_bc = top_bc / (-1.0 * Physics_Base.q_e); this.bottom_bc = bottom_bc / (-1.0 * Physics_Base.q_e);

            // and get the corresponding permittivities
            top_eps = Geom_Tool.GetLayer(layers, top_pos).Permitivity;
            bottom_eps = Geom_Tool.GetLayer(layers, bottom_pos).Permitivity;

            if (flexpde_inputfile != null)
                Create_FlexPDE_Input_File(flexpde_inputfile, layers);
        }
    }
}