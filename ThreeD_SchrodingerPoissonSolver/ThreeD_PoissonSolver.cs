using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace ThreeD_SchrodingerPoissonSolver
{
    public class ThreeD_PoissonSolver: FlexPDE_Base
    {
        bool natural_topbc = false;

        Experiment exp;
        string dens_filename = "car_dens.dat";
        string densdopent_filename = "dens_3D_dopents.dat";
        string densderiv_filename = "rho_prime.dat";
        string pot_filename = "phi.dat";
        string new_pot_filename = "new_phi.dat";
        string xc_pot_calc_filename = "xc_pot_calc.dat";

        string gphi_filename = "gphi.dat";

        double top_bc, split_bc1, split_bc2, bottom_bc;
        double surface;

        double split_width, split_length, top_length;

        double z_2DEG;

        double y_scaling, z_scaling;
        double t = 0.0;

        public ThreeD_PoissonSolver(Experiment exp, bool using_external_code, Dictionary<string, object> input)
            : base(using_external_code, input)
        {
            this.exp = exp;
            // check if the top layer is air... if so, we need to use natural boundary conditions on the upper surface (this is almost always going to be the case in 3D)
            natural_topbc = (exp.Layers[exp.Layers.Length - 1].Material == Material.Air);

            // calculate scaling factors (x is used as the reference dimension)
            // w = a_y * y such that the y dimension has the same length as the x dimension
            y_scaling = (exp.Nx_Pot * exp.Dx_Pot) / (exp.Ny_Pot * exp.Dy_Pot);
            // w = a_z * z such that the z dimension has the same length as the x dimension
            z_scaling = (exp.Nx_Pot * exp.Dx_Pot) / (exp.Nz_Pot * exp.Dz_Pot);
        }

        protected override Band_Data Parse_Potential(string[] data)
        {
            string[] new_data = Trim_Potential_File(data);
            return Band_Data.Parse_Band_Data(new_data, exp.Nx_Dens, exp.Ny_Dens, exp.Nz_Dens);
        }

        protected override void Save_Data(Band_Data density, string input_file_name)
        {
            // check that the data isn't all zeros and if it is, add a small 1e-8 perturbation at the centre
            if (density.Min() > -1e-8 && density.Max() < 1e-8)
                density.vol[density.vol.Length / 2][(int)(density.vol[density.vol.Length / 2].Rows / 2), (int)(density.vol[density.vol.Length / 2].Cols / 2)] = 1e-8;

            density.Save_3D_Data(input_file_name, exp.Dx_Dens, exp.Dy_Dens * y_scaling, exp.Dz_Dens * z_scaling, exp.Xmin_Dens, exp.Ymin_Dens * y_scaling, exp.Zmin_Dens * z_scaling);
        }

        /// <summary>
        /// calculates the Laplacian in the plane of the 2DEG.  
        /// assumes that the charge density solves the Poisson equation in the growth direction
        /// NOTE! this is not scaled as this is not (or at least should not) be used inside ThreeD_PoissonSolver
        /// </summary>
        public override Band_Data Calculate_Laplacian(Band_Data input_data)
        {
            DoubleMatrix[] result = new DoubleMatrix[exp.Nz_Dens];
            for (int k = 0; k < exp.Nz_Dens; k++)
                result[k] = new DoubleMatrix(exp.Nx_Dens, exp.Ny_Dens);
            DoubleMatrix[] data = input_data.vol;

            for (int k = 1; k < exp.Nz_Dens - 1; k++)
                for (int i = 1; i < exp.Nx_Dens - 1; i++)
                    for (int j = 1; j < exp.Ny_Dens - 1; j++)
                    {
                        double pos_x = i * exp.Dx_Dens + exp.Xmin_Dens;
                        double pos_y = j * exp.Dy_Dens + exp.Ymin_Dens;
                        double pos_z = k * exp.Dz_Dens + exp.Zmin_Dens;

                        // the factors multiplying the Laplacian in the longitudinal direction
                        double factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_x + 0.5 * exp.Dx_Dens, pos_y, pos_z).Permitivity / (exp.Dx_Dens * exp.Dx_Dens);
                        double factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_x - 0.5 * exp.Dx_Dens, pos_y, pos_z).Permitivity / (exp.Dx_Dens * exp.Dx_Dens);
                        result[k][i, j] = (factor_minus * data[k][i - 1, j] + factor_plus * data[k][i + 1, j] - (factor_plus + factor_minus) * data[k][i, j]);

                        // and in the transverse direction
                        factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_x, pos_y + 0.5 * exp.Dy_Dens, pos_z).Permitivity / (exp.Dy_Dens * exp.Dy_Dens);
                        factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_x, pos_y - 0.5 * exp.Dy_Dens, pos_z).Permitivity / (exp.Dy_Dens * exp.Dy_Dens);
                        result[k][i, j] += (factor_minus * data[k][i, j - 1] + factor_plus * data[k][i, j + 1] - (factor_plus + factor_minus) * data[k][i, j]);

                        // and in the growth direction
                        factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_x, pos_y, pos_z + 0.5 * exp.Dz_Dens).Permitivity / (exp.Dz_Dens * exp.Dz_Dens);
                        factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_x, pos_y, pos_z - 0.5 * exp.Dz_Dens).Permitivity / (exp.Dz_Dens * exp.Dz_Dens);
                        result[k][i, j] += (factor_minus * data[k - 1][i, j] + factor_plus * data[k + 1][i, j] - (factor_plus + factor_minus) * data[k][i, j]);
                    }

            return new Band_Data(result);
        }

        public override void Initiate_Poisson_Solver(Dictionary<string, double> device_dimension, Dictionary<string, double> boundary_conditions)
        {
            if (natural_topbc)
                boundary_conditions["top_V"] = 0.0;

            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // with a factor to convert from V to meV zC^-1
            top_bc = boundary_conditions["top_V"] * Physics_Base.energy_V_to_meVpzC;
            split_bc1 = boundary_conditions["split_V1"] * Physics_Base.energy_V_to_meVpzC;
            split_bc2 = boundary_conditions["split_V2"] * Physics_Base.energy_V_to_meVpzC;
            bottom_bc = boundary_conditions["bottom_V"] * Physics_Base.energy_V_to_meVpzC;

            // and save the split width and surface charge
            this.top_length = device_dimension["top_length"];
            this.split_width = device_dimension["split_width"];
            this.split_length = device_dimension["split_length"];
            this.surface = boundary_conditions["surface"];

            // and find the interface depth (for plotting)
            this.z_2DEG = device_dimension["interface_depth"] - 5.0;

            if (flexpde_script != null)
                Create_FlexPDE_File(top_bc, top_length, split_bc1, split_bc2, split_width, split_length, surface, bottom_bc, flexpde_script);
        }

        protected override Band_Data Get_ChemPot_From_External(Band_Data density)
        {
            Save_Data(density, dens_filename);

            return Get_Data_From_External(initcalc_location, flexpde_options, initcalc_result_filename);
        }

        public void Create_FlexPDE_File(double top_bc, double top_length, double split_bc1, double split_bc2, double split_width, double split_length, double surface, double bottom_bc, string output_file)
        {
            StreamWriter sw = new StreamWriter(output_file);

            // write out output file
            sw.WriteLine("TITLE \'Full Split Gate Geometry\'");
            sw.WriteLine("COORDINATES cartesian3");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            sw.WriteLine("\tERRLIM=1e-2");
            sw.WriteLine("\tGRIDLIMIT=20");
            sw.WriteLine("DEFINITIONS");
            sw.WriteLine("\trho");
            sw.WriteLine("\tband_gap");
            sw.WriteLine();
            // and the tables for carrier and donor densities
            sw.WriteLine("\trho_carrier = TABLE(\'" + dens_filename + "\', x, y, z)");
            sw.WriteLine("\trho_dopent = TABLE(\'" + densdopent_filename + "\', x, y, z)");
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\tlx = " + (exp.Dx_Pot * exp.Nx_Pot).ToString());
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Scale factors");
            sw.WriteLine("\ty_scaling = " + y_scaling.ToString());
            sw.WriteLine("\tz_scaling = " + z_scaling.ToString());
            sw.WriteLine();
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = " + surface.ToString() + " * z_scaling");
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in meV zC^-1)");
            sw.WriteLine("\tsplit_V1 = " + split_bc1.ToString());
            sw.WriteLine("\tsplit_V2 = " + split_bc2.ToString());
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString() + "");
            sw.WriteLine("\tsplit_length = " + split_length.ToString());
            sw.WriteLine("\ttop_length = " + top_length.ToString());
            sw.WriteLine("\tsplit_depth = 10! depth of the split gate metal material");
            sw.WriteLine("\ttop_depth = 10! depth of the top gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + z_2DEG.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivity");
            sw.WriteLine("\teps");
            sw.WriteLine();
            // other physical parameters
            sw.WriteLine("\tq_e = " + Physics_Base.q_e.ToString() + "! charge of electron in zC");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: dx(eps * dx(u)) + y_scaling * dy(eps * y_scaling * dy(u)) + z_scaling * dz(eps * z_scaling * dz(u)) = - rho\t! Poisson's equation");
            sw.WriteLine();
            sw.WriteLine("EXTRUSION");
            sw.WriteLine("\tSURFACE \"Substrate\"\tz = " + exp.Layers[0].Zmax.ToString() + " * z_scaling");
            int layercount = 1;
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                sw.WriteLine("\tSURFACE	\"" + layercount.ToString() + "\"\tz = " + exp.Layers[i].Zmax.ToString() + " * z_scaling");
                layercount++;

                if (exp.Layers[i].Layer_No == Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No)
                {
                    sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                    sw.WriteLine("\tSURFACE \"" + layercount.ToString() + "\"\tz = split_depth * z_scaling");
                    layercount++;
                }
                else if (exp.Layers[i].Material == Material.PMMA)
                {
                    sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                    sw.WriteLine("\tSURFACE \"" + layercount.ToString() + "\"\tz = " + exp.Layers[i].Zmax.ToString() + " * z_scaling + top_depth * z_scaling");
                    layercount++;
                }
            }
            sw.WriteLine();
            sw.WriteLine("BOUNDARIES");
            sw.WriteLine("\tSURFACE \"Substrate\"	VALUE(u) = bottom_bc");
            sw.WriteLine("\tSURFACE \"" + (Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No - 1).ToString() + "\" NATURAL(u) = surface_bc");
            sw.WriteLine("\tSURFACE \"" + (layercount - 1).ToString() + "\" NATURAL(u) = 0");
            sw.WriteLine();
            sw.WriteLine("\tREGION 1");
            layercount = 1;
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                layercount++;

                if (exp.Layers[i].Layer_No <= Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No)
                    sw.WriteLine("\t\trho = rho_carrier + rho_dopent");
                else
                    sw.WriteLine("\t\trho = 0.0");

                sw.WriteLine("\t\teps = " + exp.Layers[i].Permitivity.ToString());
                sw.WriteLine("\t\tband_gap = " + exp.Layers[i].Band_Gap.ToString());
                sw.WriteLine();

                if (exp.Layers[i].Layer_No == Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No || exp.Layers[i].Material == Material.PMMA)
                {
                    sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                    sw.WriteLine("\t\trho = 0.0");
                    sw.WriteLine("\t\teps = " + exp.Layers[i + 1].Permitivity.ToString());
                    sw.WriteLine("\t\tband_gap = " + exp.Layers[i + 1].Band_Gap.ToString());
                    sw.WriteLine();
                    layercount++;
                }
            }
            int max_layers = layercount - 1;
            sw.WriteLine("\t\tSTART(-lx / 2, -ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (-lx / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (lx / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (lx / 2, -ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 2 ! left split gate");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V1");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V1");
            sw.WriteLine("\t\tLAYER \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-split_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (split_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tVALUE(u) = split_V1");
            sw.WriteLine("\t\tLINE TO (split_length / 2, split_width / 2 * y_scaling) TO (-split_length / 2, split_width / 2 * y_scaling) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 3 ! right split gate");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V2");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V2");
            sw.WriteLine("\t\tLAYER \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-split_length / 2, -ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (split_length / 2, -ly / 2 * y_scaling)");
            sw.WriteLine("\t\tVALUE(u) = split_V2");
            sw.WriteLine("\t\tLINE TO (split_length / 2, -split_width / 2 * y_scaling) TO (-split_length / 2, -split_width / 2 * y_scaling) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 4 ! top gate");
            sw.WriteLine("\t\tSURFACE \"" + (max_layers - 2).ToString() + "\" VALUE(u) = top_V");
            sw.WriteLine("\t\tSURFACE \"" + (max_layers - 1).ToString() + "\" VALUE(u) = top_V");
            sw.WriteLine("\t\tLAYER \"" + (max_layers - 1).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-top_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (top_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tVALUE(u) = top_V");
            sw.WriteLine("\t\tLINE TO (top_length / 2, -ly / 2 * y_scaling) TO (-top_length / 2, -ly / 2 * y_scaling) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\t\tFRONT(z - well_depth * z_scaling, 20 * z_scaling)");
            sw.WriteLine();
            sw.WriteLine("!MONITORS");
            sw.WriteLine("\t!CONTOUR(rho) ON z = well_depth * z_scaling");
            sw.WriteLine("\t!CONTOUR(u) ON z = well_depth * z_scaling");
            sw.WriteLine("\t!CONTOUR(u) ON x = 0");
            sw.WriteLine("\t!CONTOUR(u) ON y = 0");
            sw.WriteLine("PLOTS");
            sw.WriteLine("\t!CONTOUR(rho) ON x = 0");
            sw.WriteLine("\t!CONTOUR(u) ON x = 0");
            sw.WriteLine("\t!CONTOUR(u) ON y = 0");
            sw.WriteLine("\t!CONTOUR(rho) ON z = well_depth * z_scaling");
            sw.WriteLine("\t!CONTOUR(- q_e * u + 0.5 * band_gap) ON z = well_depth * z_scaling");
            sw.WriteLine("\t!ELEVATION(rho) FROM (0,0, " + Geom_Tool.Get_Zmin(exp.Layers).ToString() + " * z_scaling) TO (0, 0, " + exp.Layers[exp.Layers.Length - 1].Zmax.ToString() + " * z_scaling)");
            sw.WriteLine("\t!ELEVATION(- q_e * u + 0.5 * band_gap) FROM (0, 0, " + Geom_Tool.Get_Zmin(exp.Layers).ToString() + " * z_scaling) TO (0, 0, " + exp.Layers[exp.Layers.Length - 1].Zmax.ToString() + " * z_scaling)");
            sw.WriteLine("\t!ELEVATION(- q_e * u + 0.5 * band_gap) FROM (0, -ly / 2 * y_scaling, well_depth * z_scaling) TO (0, ly / 2 * y_scaling, well_depth * z_scaling)");
            sw.WriteLine("\t!ELEVATION(- q_e * u + 0.5 * band_gap) FROM (-lx/2, 0, well_depth * z_scaling) TO (lx / 2, 0, well_depth * z_scaling)");
            sw.WriteLine();
            sw.WriteLine("\tTABLE(u) ZOOM (" + exp.Xmin_Dens.ToString() + ", " + (y_scaling * exp.Ymin_Dens).ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Nx_Dens - 1) * exp.Dx_Dens).ToString() + ", " + (y_scaling * (exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Nx_Dens.ToString() + ", " + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"pot.dat\"");
            sw.WriteLine("\tTRANSFER(u) FILE = \'" + pot_filename + "\'");
            sw.WriteLine("\tTRANSFER(0.0 * u) FILE = \'" + new_pot_filename + "\' ! dummy file for smoother function");
            sw.WriteLine();
            sw.WriteLine("END");

            sw.Close();
        }
        
        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi)
        {
            Save_Data(gphi, gphi_filename);
            Save_Data(rho_prime.Spin_Summed_Data, densderiv_filename);
            Create_NewtonStep_File(top_length, split_width, split_length, flexpde_script, T);

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

        public void Create_NewtonStep_File(double top_length, double split_width, double split_length, string output_file, double t)
        {
            StreamWriter sw = new StreamWriter(output_file);

            // generate summary containing residual integrals for car dens (will be used later for x and phi)
            string window_function_string;
            if (exp.Zmin_Dens < 0.0)
                window_function_string = "(ustep(z + " + (-1.0 * z_scaling * exp.Zmin_Dens).ToString() + ")";
            else
                window_function_string = "(ustep(z - " + (z_scaling * exp.Zmin_Dens).ToString() + ")";
            double zmax = z_scaling * (exp.Zmin_Dens + exp.Nz_Dens * exp.Dz_Dens);
            if (zmax < 0.0)
                window_function_string = window_function_string + " - ustep(z + " + (-1.0 * zmax).ToString() + "))";
            else
                window_function_string = window_function_string + " - ustep(z - " + zmax.ToString() + "))";

            string minus_g_phi = "-1.0 * gphi * " + window_function_string;

            // write out output file
            sw.WriteLine("TITLE \'Full Split Gate Geometry\'");
            sw.WriteLine("COORDINATES cartesian3");
            sw.WriteLine();
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine();
            sw.WriteLine("SELECT");
            sw.WriteLine("\tERRLIM=1e-2");
            sw.WriteLine("\tGRIDLIMIT=20");
            sw.WriteLine();
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\tTRANSFERMESH(\'" + pot_filename + "\', phi)");
            sw.WriteLine("\tTRANSFER(\'" + new_pot_filename + "\', new_phi)");
            sw.WriteLine();
            // and the tables for carrier and donor densities
            //sw.WriteLine("\tcar_dens = SMOOTH(" + exp.Dy_Dens.ToString() + ") TABLE(\'" + dens_filename + "\')");
            sw.WriteLine("\tgphi = TABLE(\'" + gphi_filename + "\')");
            sw.WriteLine("\txc_pot = TABLE(\'" + xc_pot_filename + "\')");
            sw.WriteLine("\txc_pot_calc = TABLE(\'xc_pot_calc.dat\')");
            sw.WriteLine("\tcar_dens = TABLE(\'" + dens_filename + "\')");
            sw.WriteLine("\trho_prime = TABLE(\'" + densderiv_filename + "\')");
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\ty_scaling = " + y_scaling.ToString());
            sw.WriteLine("\tz_scaling = " + z_scaling.ToString());
            sw.WriteLine("\tlx = " + (exp.Dx_Pot * exp.Nx_Pot).ToString());
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine("\tlz = " + (exp.Dz_Pot * exp.Nz_Pot).ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = 0.0");
            sw.WriteLine("\tsurface_bc = 0.0");
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in meV zC^-1)");
            sw.WriteLine("\tsplit_V = 0.0");
            sw.WriteLine("\ttop_V = 0.0");
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString() + "");
            sw.WriteLine("\tsplit_length = " + split_length.ToString());
            sw.WriteLine("\ttop_length = " + top_length.ToString());
            sw.WriteLine("\tsplit_depth = 10! depth of the split gate metal material");
            sw.WriteLine("\ttop_depth = 10! depth of the top gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + z_2DEG.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivity");
            sw.WriteLine("\teps");
            sw.WriteLine();
            // other physical parameters
            sw.WriteLine("\tq_e = " + Physics_Base.q_e.ToString() + "! charge of electron in zC");
            sw.WriteLine();
            sw.WriteLine("\tt = " + t.ToString() + "! Mixing parameter from previous iteration");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: -1.0 * dx(eps * dx(u)) - 1.0 * y_scaling * dy(eps * y_scaling * dy(u)) - 1.0 * z_scaling * dz(eps * z_scaling * dz(u)) - rho_prime * u = " + minus_g_phi + " \t! Poisson's equation");
            sw.WriteLine();
            sw.WriteLine("EXTRUSION");
            sw.WriteLine("\tSURFACE \"Substrate\"\tz = " + exp.Layers[0].Zmax.ToString() + " * z_scaling");
            int layercount = 1;
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                sw.WriteLine("\tSURFACE	\"" + layercount.ToString() + "\"\tz = " + exp.Layers[i].Zmax.ToString() + " * z_scaling");
                layercount++;

                if (exp.Layers[i].Layer_No == Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No)
                {
                    sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                    sw.WriteLine("\tSURFACE \"" + layercount.ToString() + "\"\tz = split_depth * z_scaling");
                    layercount++;
                }
                else if (exp.Layers[i].Material == Material.PMMA)
                {
                    sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                    sw.WriteLine("\tSURFACE \"" + layercount.ToString() + "\"\tz = " + exp.Layers[i].Zmax.ToString() + " * z_scaling + top_depth * z_scaling");
                    layercount++;
                }
            }
            sw.WriteLine();
            sw.WriteLine("BOUNDARIES");
            sw.WriteLine("\tSURFACE \"Substrate\"	VALUE(u) = 0");
            sw.WriteLine("\tSURFACE \"" + (layercount - 1).ToString() + "\" NATURAL(u) = 0");
            sw.WriteLine();
            sw.WriteLine("\tREGION 1");
            layercount = 1;
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                layercount++;

                sw.WriteLine("\t\teps = " + exp.Layers[i].Permitivity.ToString());
                sw.WriteLine();

                if (exp.Layers[i].Layer_No == Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No || exp.Layers[i].Material == Material.PMMA)
                {
                    sw.WriteLine("\t\tLAYER \"" + layercount.ToString() + "\"");
                    sw.WriteLine("\t\teps = " + exp.Layers[i + 1].Permitivity.ToString());
                    sw.WriteLine();
                    layercount++;
                }
            }
            int max_layers = layercount - 1;
            sw.WriteLine("\t\tSTART(-lx / 2, -ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (-lx / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (lx / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (lx / 2, -ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 2 ! left split gate");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tLAYER \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-split_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (split_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (split_length / 2, split_width / 2 * y_scaling) TO (-split_length / 2, split_width / 2 * y_scaling) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 3 ! right split gate");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tSURFACE \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tLAYER \"" + (Geom_Tool.Find_Layer_Above_Surface(exp.Layers).Layer_No - 1).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-split_length / 2, -ly / 2 *  y_scaling)");
            sw.WriteLine("\t\tLINE TO (split_length / 2, -ly / 2 *  y_scaling)");
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (split_length / 2, -split_width / 2 *  y_scaling) TO (-split_length / 2, -split_width / 2 *  y_scaling) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 4 ! top gate");
            sw.WriteLine("\t\tSURFACE \"" + (max_layers - 2).ToString() + "\" VALUE(u) = top_V");
            sw.WriteLine("\t\tSURFACE \"" + (max_layers - 1).ToString() + "\" VALUE(u) = top_V");
            sw.WriteLine("\t\tLAYER \"" + (max_layers - 1).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-top_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tLINE TO (top_length / 2, ly / 2 * y_scaling)");
            sw.WriteLine("\t\tVALUE(u) = top_V");
            sw.WriteLine("\t\tLINE TO (top_length / 2, -ly / 2 * y_scaling) TO (-top_length / 2, -ly / 2 * y_scaling) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tRESOLVE(" + minus_g_phi + ")");
            sw.WriteLine();
            sw.WriteLine("!MONITORS");
            sw.WriteLine("\t!CONTOUR(rho) ON z = well_depth * z_scaling");
            sw.WriteLine("\t!CONTOUR(u) ON z = well_depth * z_scaling");
            sw.WriteLine("\t!CONTOUR(u) ON x = 0");
            sw.WriteLine("\t!CONTOUR(u) ON y = 0");
            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(car_dens) ON z = well_depth * z_scaling ON GRID(x, y / y_scaling)");
            sw.WriteLine("\tCONTOUR(1424.0 * 0.5 - q_e * (phi + t * new_phi)) ON z = well_depth * z_scaling ON GRID(x, y / y_scaling)");
            sw.WriteLine("\tCONTOUR(u * q_e) ON z = well_depth * z_scaling ON GRID(x, y / y_scaling)");
            sw.WriteLine("\tCONTOUR(car_dens) ON x = 0 ON GRID(y / y_scaling, z / z_scaling)");
            sw.WriteLine("\tCONTOUR(u * q_e) ON x = 0 ON GRID(y / y_scaling, z / z_scaling)");
            sw.WriteLine("\tCONTOUR(u * q_e) ON y = 0 ON GRID(x, z / z_scaling)");
            sw.WriteLine("\tCONTOUR(" + minus_g_phi + ") ON x = 0 ON GRID(y / y_scaling, z / z_scaling)");
            sw.WriteLine("\tCONTOUR(" + minus_g_phi + ") ON z = well_depth * z_scaling ON GRID(x, y / y_scaling)");
            sw.WriteLine("\tCONTOUR(rho_prime) ON z = well_depth * z_scaling ON GRID(x, y / y_scaling)");
            sw.WriteLine();
            sw.WriteLine("\tTRANSFER(phi + t * new_phi) FILE = \'" + pot_filename + "\'");
            sw.WriteLine("\tTRANSFER(u) FILE = \'" + new_pot_filename + "\'");
            sw.WriteLine();
            sw.WriteLine("\tTABLE(phi + t * new_phi) ZOOM (" + exp.Xmin_Dens.ToString() + ", " + (y_scaling * exp.Ymin_Dens).ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Nx_Dens - 1) * exp.Dx_Dens).ToString() + ", " + (y_scaling * (exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Nx_Dens.ToString() + ", " + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"y.dat\"");
            sw.WriteLine("\tTABLE(u) ZOOM (" + exp.Xmin_Dens.ToString() + ", " + (y_scaling * exp.Ymin_Dens).ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Nx_Dens - 1) * exp.Dx_Dens).ToString() + ", " + (y_scaling * (exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Nx_Dens.ToString() + ", " + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"x.dat\"");
            sw.WriteLine();
            sw.WriteLine("END");

            sw.Close();
        }

        public override Band_Data Chemical_Potential
        {
            get
            {
                string[] lines = File.ReadAllLines("y.dat");
                string[] data = Trim_Potential_File(lines);

                // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
                return Physics_Base.q_e * Parse_Potential(data);
            }
        }
    }
}
