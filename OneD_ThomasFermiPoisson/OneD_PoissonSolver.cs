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
        Experiment exp;

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

        public OneD_PoissonSolver(Experiment exp, bool using_flexPDE, string flexPDE_input, string flexPDE_location, double tol)
            : base(using_flexPDE, flexPDE_input, flexPDE_location, tol)
        {
            this.exp = exp;

            // generate Laplacian matrix (spin-resolved)
            if (!flexPDE)
            {
                // first check whether the grids for the potential and the density are the same
                //if (exp.Nz_Dens != exp.Nz_Pot || exp.Dz_Dens != exp.Dz_Pot || exp.Zmin_Dens != exp.Zmin_Pot)
                //    throw new Exception("Error - when not using flexPDE, the density and potential grids must be the same!");

                laplacian = Generate_Laplacian(exp.Layers);
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

        protected override Band_Data Parse_Potential(string[] data)
        {
            // and check that there is the right number of data points back
            if (data.Length != exp.Nz_Pot)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleVector(exp.Nz_Pot));
            for (int i = 0; i < exp.Nz_Pot; i++)
                result.vec[i] = double.Parse(data[i]);

            return result;
        }

        protected override Band_Data Get_BandEnergy_On_Regular_Grid(Band_Data charge_density)
        {
            // set the top and bottom boundary conditions where [0] is the bottom of the device
            charge_density.vec[0] = bottom_bc * -1.0 * bottom_eps / (exp.Dz_Pot * exp.Dz_Pot);
            charge_density.vec[charge_density.Length - 1] = top_bc * -1.0 * top_eps / (exp.Dz_Pot * exp.Dz_Pot);

            // solve Poisson's equation
            Band_Data potential = new Band_Data(lu_fact.Solve(charge_density.vec));

            // return band energy
            return -1.0 * Physics_Base.q_e * potential * 6.2415093;  // the factor of 6.24 is because 1 zC V = 6.24 meV
        }

        /// <summary>
        /// Generates a laplacian matrix in one-dimension on a regular grid with Dirichlet BCs and assuming constant permitivity
        /// </summary>
        DoubleMatrix Generate_Laplacian()
        {
            // the factor which multiplies the Laplace equation
            double factor = -1.0 * Physics_Base.epsilon / (exp.Dz_Pot * exp.Dz_Pot);

            DoubleMatrix result = new DoubleMatrix(exp.Nz_Pot, exp.Nz_Pot);
            for (int i = 0; i < exp.Nz_Pot - 1; i++)
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
            result[exp.Nz_Pot - 1, exp.Nz_Pot - 1] = 1.0 * factor;
            result[exp.Nz_Pot - 1, exp.Nz_Pot - 2] = 0.0;

            return result;
        }

        /// <summary>
        /// Generates a laplacian matrix in one-dimension on a regular grid with Dirichlet BCs wot varying permitivity
        /// </summary>
        DoubleMatrix Generate_Laplacian(ILayer[] layers)
        {
            DoubleMatrix result = new DoubleMatrix(exp.Nz_Pot, exp.Nz_Pot);
            double factor_plus, factor_minus;

            // cycle through the structure and fill in the Laplacian with the correct permittivities
            for (int i = 1; i < exp.Nz_Pot - 1; i++)
            {
                double eps_plus = Geom_Tool.GetLayer(layers, i * exp.Dz_Pot + 0.5 * exp.Dz_Pot + exp.Zmin_Pot).Permitivity;
                double eps_minus = Geom_Tool.GetLayer(layers, i * exp.Dz_Pot - 0.5 * exp.Dz_Pot + exp.Zmin_Pot).Permitivity;

                // the factor which multiplies the Laplace equation
                factor_plus = -1.0 * eps_plus / (exp.Dz_Pot * exp.Dz_Pot);
                factor_minus = -1.0 * eps_minus / (exp.Dz_Pot * exp.Dz_Pot);

                // on-diagonal term
                result[i, i] = factor_minus + factor_plus;
                // off-diagonal
                result[i, i - 1] = -1.0 * factor_minus;
                result[i, i + 1] = -1.0 * factor_plus;
            }

            // and fix boundary conditions
            double factor = -1.0 * Geom_Tool.GetLayer(layers, exp.Zmin_Pot).Permitivity / (exp.Dz_Pot * exp.Dz_Pot);
            result[0, 0] = 1.0 * factor;
            result[0, 1] = 0.0;
            factor = -1.0 * Geom_Tool.GetLayer(layers, (exp.Nz_Pot - 1) * exp.Dz_Pot + exp.Zmin_Pot).Permitivity / (exp.Dz_Pot * exp.Dz_Pot);
            result[exp.Nz_Pot - 1, exp.Nz_Pot - 1] = 1.0 * factor;
            result[exp.Nz_Pot - 1, exp.Nz_Pot - 2] = 0.0;

            return result;
        }

        public double Get_Surface_Charge(Band_Data band_offset, ILayer[] layers)
        {
            // calculate the electric field just below the surface
            int surface = (int)(-1.0 * Math.Floor(Geom_Tool.Get_Zmin(layers) / exp.Dz_Pot));
            double eps = layers[Geom_Tool.Find_Layer_Below_Surface(layers)].Permitivity;
            // by Gauss' theorem, rho = - epsilon_0 * epsilon_r * dV/dz
            double surface_charge = eps * (band_offset[surface] - band_offset[surface - 1]) / exp.Dz_Pot;
            // divide by q_e to convert the band energy into a potential
            surface_charge /= Physics_Base.q_e;
            // and divide by dz to give a density
            surface_charge /= exp.Dz_Pot;

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
            // Poisson's equation
            sw.WriteLine("\tu: div(eps * grad(u)) = -rho\t! Poisson's equation");
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

        protected void Create_FlexPDE_Input_File(string flexPDE_input)
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
            sw.WriteLine("\tband_gap");
            // number of lattice sites that the density needs to be output to
            sw.WriteLine("\tnz = " + exp.Nz_Pot.ToString());
            // size of the sample
            sw.WriteLine("\tlz = " + (exp.Nz_Pot * exp.Dz_Pot).ToString());
            // the top boundary condition on the surface of the sample
            sw.WriteLine("\t! Boundary conditions");
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine("\tbottom_V = " + bottom_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivities");
            sw.WriteLine("\teps_0 = " + Physics_Base.epsilon_0.ToString());
            // relative permitivity of materials
            sw.WriteLine("\teps_r_GaAs = " + Physics_Base.epsilon_r_GaAs.ToString());
            sw.WriteLine("\teps_r_AlGaAs = " + Physics_Base.epsilon_r_AlGaAs.ToString());
            sw.WriteLine("\teps_pmma = " + Physics_Base.epsilon_pmma.ToString());
            sw.WriteLine("\teps");
            sw.WriteLine();
            // other physical parameters
            sw.WriteLine("\tq_e = " + Physics_Base.q_e.ToString() + "! charge of electron in zC");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: div(eps * grad(u)) = -rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");

            // cycle through layers
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\tREGION " + exp.Layers[i].Layer_No.ToString());
                sw.WriteLine("\t\trho = TABLE(\'" + dens_filename + "\', x)");
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(exp.Layers[i].Material));
                sw.WriteLine("\t\tband_gap = " + exp.Layers[i].Band_Gap.ToString());
                sw.WriteLine("\t\tSTART(" + exp.Layers[i].Zmin.ToString() + ")");
                if (i == 1)
                    sw.WriteLine("\t\tPOINT VALUE(u) = top_V");
                sw.WriteLine("\t\tLINE TO (" + exp.Layers[i].Zmax.ToString() + ")");
                if (i == exp.Layers.Length - 1)
                    sw.WriteLine("\t\tPOINT VALUE(u) = bottom_V");
                sw.WriteLine();
            }

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tELEVATION(-q_e * u + 0.5 * band_gap) FROM (-lz) TO (0)");
            sw.WriteLine("\tELEVATION(rho) FROM (-lz) TO (0)");
            sw.WriteLine("\tELEVATION(u) FROM (" + exp.Zmin_Pot.ToString() + ") TO (" + (exp.Zmin_Pot + exp.Nz_Pot * exp.Dz_Pot).ToString() + ") EXPORT(nz) FORMAT \'#1\' FILE=\'pot.dat\'");
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
                Create_FlexPDE_Input_File(flexpde_inputfile);
        }

        protected override void Save_Density_Data(Band_Data density, string input_file_name)
        {
            density.Save_1D_Data(input_file_name, exp.Dz_Dens, exp.Zmin_Dens);
        }
    }
}