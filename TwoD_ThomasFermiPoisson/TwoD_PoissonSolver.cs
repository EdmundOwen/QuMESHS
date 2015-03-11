using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using CenterSpace.NMath.Core;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace TwoD_ThomasFermiPoisson
{
    public class TwoD_PoissonSolver : FlexPDE_Base
    {
        Experiment exp;
        string dens_filename = "car_dens.dat";
        string densdopent_filename = "dens_2D_dopents.dat";
        string densderiv_filename = "rho_prime.dat";
        string pot_filename = "phi.dat";
        string new_pot_filename = "new_phi.dat";
        string xc_pot_calc_filename = "xc_pot_calc.dat";

        string gphi_filename = "gphi.dat";

        string chempot_result_filename = "y.dat";

        double top_bc, split_bc, bottom_bc;
        double split_width;
        double surface;

        public TwoD_PoissonSolver(Experiment exp, bool using_external_code, Dictionary<string, object> input)
            : base(using_external_code, input)
        {
            this.exp = exp;
        }

        protected override Band_Data Parse_Potential(string[] data)
        {
            string[] new_data = Trim_Potential_File(data);
            return Band_Data.Parse_Band_Data(new_data, exp.Ny_Dens, exp.Nz_Dens);
        }

        public Band_Data Calculate_Point_Potential(double top_bc, double split_bc, double split_width, double surface, double bottom_bc, int i, int j, double U)
        {
            // don't bother if on the edge of the domain
            if (i == 0 || i == exp.Ny_Dens - 1 || j == 0 || j == exp.Nz_Dens - 1)
                return null;

            double ypos = exp.Ymin_Dens + i * exp.Dy_Dens;
            double zpos = exp.Zmin_Dens + j * exp.Dz_Dens;

            string dens_line = "\trho_carrier = " + U.ToString() + " * min( (1 - abs(x - (" + ypos.ToString() + ")) / " + exp.Dy_Dens.ToString() + ") , "
                                    + "(1 - abs(y  - (" + zpos.ToString() + ")) / " + exp.Dz_Dens.ToString() + ")) * "
                                    + "upulse(x - (" + (ypos - exp.Dy_Dens).ToString() + "), x - (" + (ypos + exp.Dy_Dens).ToString() + ")) * "
                                    + "upulse(y - (" + (zpos - exp.Dz_Dens).ToString() + "), y - (" + (zpos + exp.Dz_Dens).ToString() + "))"
                                    + "\n\trho_dopent = 0";

            // create flexPDE file for a point charge at (ypos, zpos)
            Create_FlexPDE_File(top_bc, split_bc, split_width, surface, bottom_bc, initcalc_result_filename, dens_line);

            // get the data using flexPDE
            return Get_Data_From_External(initcalc_location, flexpde_options, initcalc_result_filename);
        }

        public Band_Data Get_Chemical_Potential(Band_Data background, Band_Data[] point_potentials, Band_Data carrier_densities, double norm)
        {
            Band_Data result = new Band_Data(new DoubleMatrix(exp.Ny_Dens, exp.Nz_Dens));

            // cycle over each point and sum the potentials due to the point-wise carrier densities
            for (int i = 0; i < exp.Ny_Dens; i++)
                for (int j = 0; j < exp.Nz_Dens; j++)
                {
                    // add the background potential
                    result.mat[i, j] = background.mat[i, j];

                    for (int k = 0; k < exp.Ny_Dens * exp.Nz_Dens; k++)
                    {
                        // don't do anything if k is on the boundary (which means it's a null reference in point_potentials)
                        if (point_potentials[k] == null)
                            continue;

                        // add a weighted value of the carrier density for this point
                        result.mat[i, j] += norm * carrier_densities.mat[i, j] * point_potentials[k].mat[i, j];
                    }
                }

            // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
            return result;
        }

        protected override Band_Data Get_ChemPot_From_External(Band_Data density)
        {
            Save_Data(density, dens_filename);

            return Get_Data_From_External(initcalc_location, flexpde_options, initcalc_result_filename);
        }

        public override void Create_FlexPDE_File(double top_bc, double split_bc, double split_width, double surface, double bottom_bc, string output_file)
        {
            Create_FlexPDE_File(top_bc, split_bc, split_width, surface, bottom_bc, output_file, "rho_carrier = TABLE(\'" + dens_filename + "\', x, y)" + "\t\nrho_dopent = TABLE(\'dens_2D_dopents.dat\', x, y)");
        }

        public void Create_FlexPDE_File(double top_bc, double split_bc, double split_width, double surface, double bottom_bc, string output_file, string dens_line)
        {
            StreamWriter sw = new StreamWriter(output_file);

            sw.WriteLine("TITLE \'Split Gate\'");
            sw.WriteLine("COORDINATES cartesian2");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            sw.WriteLine("\tERRLIM=1e-5");
            sw.WriteLine("\tGRIDLIMIT=20");
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\trho");
            sw.WriteLine("\tband_gap");
            sw.WriteLine();
            // and the tables for carrier and donor densities
            sw.WriteLine(dens_line);
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine("\tlz = " + (exp.Dz_Pot * exp.Nz_Pot).ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = " + surface.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in meV zC^-1)");
            sw.WriteLine("\tsplit_V = " + split_bc.ToString());
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString());
            sw.WriteLine("\tsplit_depth = 10\t! depth of the split gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + (exp.Layers[1].Zmax - 5).ToString());
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
            sw.WriteLine("\tu: div(eps * grad(u)) = - rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");

            // cycle through layers below surface
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\tREGION " + exp.Layers[i].Layer_No.ToString());
                if (exp.Layers[i].Layer_No <= Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No)
                    sw.WriteLine("\t\trho = rho_carrier + rho_dopent");
                else
                    sw.WriteLine("\t\trho = 0.0");
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(exp.Layers[i].Material));
                sw.WriteLine("\t\tband_gap = " + exp.Layers[i].Band_Gap.ToString());

                // HACK!
                //if (i == 1)
                //    sw.WriteLine("mesh_density=0.07*exp(-0.00001*(y-well_depth)^2)*exp(-0.000001*x^2)");

                sw.WriteLine("\t\tSTART(ly / 2, " + exp.Layers[i].Zmin.ToString() + ")");
                sw.WriteLine("\t\tLINE TO (ly / 2, " + exp.Layers[i].Zmax.ToString() + ")");

                // set top gate here
                if (i == exp.Layers.Length - 1)
                    sw.WriteLine("\t\tVALUE(u) = top_V");
                // or surface condition
                if (exp.Layers[i].Zmax == 0.0)
                    sw.WriteLine("\t\tNATURAL(u) = surface_bc * upulse(x + split_width / 2 - 20, x - split_width / 2 + 20)");

                sw.WriteLine("\t\tLINE TO (-ly / 2, " + exp.Layers[i].Zmax.ToString() + ")");
                sw.WriteLine("\t\tNATURAL(u) = 0");
                sw.WriteLine("\t\tLINE TO (-ly / 2, " + exp.Layers[i].Zmin.ToString() + ")");
                // set bottom boundary condition
                if (i == 1)
                    sw.WriteLine("\t\tVALUE(u) = bottom_bc");
                sw.WriteLine("\t\tLINE TO CLOSE");
                sw.WriteLine();
            }

            // write in surface and gates
            sw.WriteLine("\tREGION " + exp.Layers.Length.ToString() + " ! Left split gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\tband_gap = 0");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(-ly / 2, 0)");
            // left split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (-ly / 2, split_depth) TO (-split_width / 2, split_depth) TO (-split_width / 2, 0) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION " + (exp.Layers.Length + 1).ToString() + "! Right split gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\tband_gap = 0");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(ly / 2, 0)");
            // right split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (ly / 2, split_depth) TO (split_width / 2, split_depth) TO (split_width / 2, 0) TO CLOSE");
            sw.WriteLine();

            // add a front commmand at the well depth to introduce a higher density of points at the interface
            sw.WriteLine("\tFRONT(y - well_depth, 1)");
            sw.WriteLine();

            // write out plotting routines
            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(- q_e * u + 0.5 * band_gap)");
            sw.WriteLine("\tCONTOUR(u)");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(u) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(rho) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(- q_e * u + 0.5 * band_gap)");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (0, " + (exp.Zmin_Dens + (exp.Nz_Dens + 1) * exp.Dz_Dens).ToString() + ") TO (0, " + (exp.Zmin_Dens - 2.0 * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (-" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth) TO (" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth)");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tELEVATION(rho) FROM (0, " + (exp.Zmin_Dens + (exp.Nz_Dens + 1) * exp.Dz_Dens).ToString() + ") TO (0, " + (exp.Zmin_Dens - 2.0 * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tELEVATION(-1.0e21*rho/q_e) FROM (-" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth) TO (" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth)");

            // and transfer the data to a file for reloading and replotting later
            sw.WriteLine();
            sw.WriteLine("\tTABLE(u) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + ((exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + initcalc_result_filename + "\"");
            sw.WriteLine("\tTRANSFER(u) FILE = \'" + pot_filename + "\'");
            sw.WriteLine("\tTRANSFER(0.0 * u) FILE = \'" + new_pot_filename  + "\' ! dummy file for smoother function");
            sw.WriteLine();

            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }

        public override void Initiate_Poisson_Solver(Dictionary<string, double> device_dimension, Dictionary<string, double> boundary_conditions)
        {
            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // with a factor to convert from V to meV zC^-1
            top_bc = boundary_conditions["top_V"] * Physics_Base.energy_V_to_meVpzC;
            split_bc = boundary_conditions["split_V"] * Physics_Base.energy_V_to_meVpzC;
            bottom_bc = boundary_conditions["bottom_V"] * Physics_Base.energy_V_to_meVpzC;

            // and save the split width and surface charge
            this.split_width = device_dimension["split_width"];
            this.surface = boundary_conditions["surface"];

            if (flexpde_script != null)
                Create_FlexPDE_File(top_bc, split_bc, split_width, surface, bottom_bc, flexpde_script);
        }

        protected override Band_Data Get_ChemPot_On_Regular_Grid(Band_Data density)
        {
            throw new NotImplementedException();
        }

        protected override void Save_Data(Band_Data density, string input_file_name)
        {
            density.Save_2D_Data(input_file_name, exp.Dy_Dens, exp.Dz_Dens, exp.Ymin_Dens, exp.Zmin_Dens);
        }


        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi)
        {
            throw new NotImplementedException();
            Save_Data(gphi, gphi_filename);
            Save_Data(rho_prime.Spin_Summed_Data, densderiv_filename);
         //   Create_NewtonStep_File(top_bc, split_bc, split_width, surface, bottom_bc, newton_filename, T);

            return Get_Data_From_External(newton_location, flexpde_options, newton_result_filename);
        }

        public Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data car_dens)
        {
            Save_Data(car_dens.Spin_Summed_Data, dens_filename);
            Save_Data(Physics_Base.Get_XC_Potential(car_dens), xc_pot_filename);

            return Calculate_Newton_Step(rho_prime, gphi);
        }

        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data car_dens, Band_Data dft_diff)
        {
            Save_Data(dft_diff + Physics_Base.Get_XC_Potential(car_dens), xc_pot_calc_filename);

            return Calculate_Newton_Step(rho_prime, gphi, car_dens);
        }

        public override void Create_NewtonStep_File(double top_bc, double split_bc, double split_width, double surface, double bottom_bc, string output_file, double t)
        {
            throw new NotImplementedException();
            StreamWriter sw = new StreamWriter(output_file);

            sw.WriteLine("TITLE \'Split Gate\'");
            sw.WriteLine("COORDINATES cartesian2");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            sw.WriteLine("\tERRLIM=1e-5");
            sw.WriteLine("\tGRIDLIMIT=20");
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\tTRANSFER(\'" + gphi_filename + "\', g_phi)");
            sw.WriteLine("\tTRANSFER(\'" + pot_filename + "\', phi)");
            sw.WriteLine();
            // and the tables for carrier and donor densities
            sw.WriteLine("\tcar_dens = TABLE(\'" + dens_filename + "\')");
            sw.WriteLine("\tdop_dens = TABLE(\'" + densdopent_filename + "\')");
            sw.WriteLine("\trho_prime = TABLE(\'" + densderiv_filename + "\')");
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine("\tlz = " + (exp.Dz_Pot * exp.Nz_Pot).ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = 0.0");// + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = 0.0");// + surface.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in meV zC^-1)");
            sw.WriteLine("\tsplit_V = 0.0");// + split_bc.ToString());
            sw.WriteLine("\ttop_V = 0.0");// + top_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString());
            sw.WriteLine("\tsplit_depth = 10\t! depth of the split gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + (exp.Layers[1].Zmax - 5).ToString());
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
            sw.WriteLine("\tu: -1.0 * div(eps * grad(u)) - rho_prime * u = -1.0 * g_phi \t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");

            // cycle through layers below surface
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\tREGION " + (exp.Layers[i].Layer_No - 1).ToString());
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(exp.Layers[i].Material));

                // HACK!
                //if (i == 1)
                //    sw.WriteLine("mesh_density=0.07*exp(-0.00001*(y-well_depth)^2)*exp(-0.000001*x^2)");

                sw.WriteLine("\t\tSTART(ly / 2, " + exp.Layers[i].Zmin.ToString() + ")");
                sw.WriteLine("\t\tLINE TO (ly / 2, " + exp.Layers[i].Zmax.ToString() + ")");

                // set top gate here
                if (i == exp.Layers.Length - 1)
                    sw.WriteLine("\t\tVALUE(u) = top_V");
                // or surface condition
                if (exp.Layers[i].Zmax == 0.0)
                    sw.WriteLine("\t\tNATURAL(u) = surface_bc * upulse(x + split_width / 2 - 20, x - split_width / 2 + 20)");

                sw.WriteLine("\t\tLINE TO (-ly / 2, " + exp.Layers[i].Zmax.ToString() + ")");
                sw.WriteLine("\t\tNATURAL(u) = 0");
                sw.WriteLine("\t\tLINE TO (-ly / 2, " + exp.Layers[i].Zmin.ToString() + ")");
                // set bottom boundary condition
                if (i == 1)
                    sw.WriteLine("\t\tVALUE(u) = bottom_bc");
                sw.WriteLine("\t\tLINE TO CLOSE");
                sw.WriteLine();
            }

            // write in surface and gates
            sw.WriteLine("\tREGION " + exp.Layers.Length.ToString() + " ! Left split gate");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(-ly / 2, 0)");
            // left split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (-ly / 2, split_depth) TO (-split_width / 2, split_depth) TO (-split_width / 2, 0) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION " + (exp.Layers.Length + 1).ToString() + "! Right split gate");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(ly / 2, 0)");
            // right split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (ly / 2, split_depth) TO (split_width / 2, split_depth) TO (split_width / 2, 0) TO CLOSE");
            sw.WriteLine();

            // write out plotting routines

            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(g_phi)");
            sw.WriteLine("\tCONTOUR(u)");
            sw.WriteLine("\tELEVATION(u) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(g_phi) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(phi)");
            sw.WriteLine("\tCONTOUR(phi * q_e) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + ((exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(car_dens) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + ((exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(car_dens + dop_dens)");
            sw.WriteLine("\tCONTOUR(u)");
            sw.WriteLine("\tCONTOUR(u * q_e) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + ((exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(g_phi) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + ((exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(g_phi)");

            // generate summary containing residual integrals for x and phi
            string window_function_string;
            if (exp.Ymin_Dens < 0.0)
                window_function_string = "(1.0 - (ustep(x + " + (-1.0 * exp.Ymin_Dens).ToString() + ")";
            else
                window_function_string= "(1.0 - (ustep(x - " + exp.Ymin_Dens.ToString() + ")";
            double xmax = exp.Ymin_Dens + exp.Ny_Dens * exp.Dy_Dens;
            if (xmax < 0.0)
                window_function_string = window_function_string + " - ustep(x + " + (-1.0 * xmax).ToString() + "))";
            else
                window_function_string = window_function_string + " - ustep(x - " + xmax.ToString() + "))";
            if (exp.Zmin_Dens < 0.0)
                window_function_string = window_function_string + " * (ustep(y + " + (-1.0 * exp.Zmin_Dens).ToString() + ")";
            else
                window_function_string = window_function_string + " * (ustep(y - " + exp.Zmin_Dens.ToString() + ")";
            double ymax = exp.Zmin_Dens + exp.Nz_Dens * exp.Dz_Dens;
            if (ymax < 0.0)
                window_function_string = window_function_string + " - ustep(y + " + (-1.0 * ymax).ToString() + ")))";
            else
                window_function_string = window_function_string + " - ustep(y - " + ymax.ToString() + "))";

            window_function_string = window_function_string + " - (ustep(y + 3) - ustep(y -13)) * (ustep(x + 353) - ustep(x + 347) - ustep(x - 353) + ustep(x - 347))";

            string residual_report_filename = "residual_g.dat";
            sw.WriteLine("\tSUMMARY EXPORT FILE = \'" + residual_report_filename + "\'");
            sw.WriteLine("\t\tREPORT ( INTEGRAL (g_phi * u * " + window_function_string + ")) AS \"residual_g_phi\"");
            sw.WriteLine("\t\tREPORT ( INTEGRAL ((div (eps * grad(u))) * u * " + window_function_string + ")) AS \"residual_g_x\"");

            // and transfer the data to a file for reloading and replotting later
            sw.WriteLine();
            sw.WriteLine("\tTRANSFER(u) FILE = \'" + new_pot_filename + "\'");
            sw.WriteLine("\tTABLE(u) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + ((exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + newton_result_filename + "\"");
            sw.WriteLine();

            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }

        public override Band_Data Calculate_Laplacian(Band_Data input_data)
        {
            DoubleMatrix result = new DoubleMatrix(exp.Ny_Dens, exp.Nz_Dens);
            DoubleMatrix data = input_data.mat;

            for (int i = 1; i < exp.Ny_Dens - 1;i++)
                for (int j = 1; j < exp.Nz_Dens - 1; j++)
                {
                    double pos_y = i * exp.Dy_Dens + exp.Ymin_Dens;
                    double pos_z = j * exp.Dz_Dens + exp.Zmin_Dens;

                    // the factors multiplying the Laplacian in the transverse direction
                    double factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_y + exp.Dy_Dens, pos_z).Permitivity / (exp.Dy_Dens * exp.Dy_Dens);
                    double factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_y - exp.Dy_Dens, pos_z).Permitivity / (exp.Dy_Dens * exp.Dy_Dens);
                    result[i, j] = (factor_minus * data[i - 1, j] + factor_plus * data[i + 1, j] - (factor_plus + factor_minus) * data[i, j]);

                    // and in the growth direction
                    factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_y, pos_z + exp.Dz_Dens).Permitivity / (exp.Dz_Dens * exp.Dz_Dens);
                    factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_y, pos_z - exp.Dz_Dens).Permitivity / (exp.Dz_Dens * exp.Dz_Dens);
                    result[i, j] += (factor_minus * data[i, j - 1] + factor_plus * data[i, j + 1] - (factor_plus + factor_minus) * data[i, j]);
                }

            return new Band_Data(result);
        }

        public override Band_Data Chemical_Potential
        { get { throw new NotImplementedException(); } }

        bool with_smoothing = false;
        public bool With_Smoothing { set { with_smoothing = value; } }
    }
}
