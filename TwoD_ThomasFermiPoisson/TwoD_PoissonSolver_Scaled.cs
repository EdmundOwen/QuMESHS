/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_PoissonSolver_Scaled : FlexPDE_Base
    {
        bool natural_topbc = false;

        Experiment exp;
        string dens_filename = "car_dens.dat";
        string densdopent_filename = "dens_2D_dopents.dat";
        string densderiv_filename = "rho_prime.dat";
        string pot_filename = "phi.dat";
        string new_pot_filename = "new_phi.dat";
        string xc_pot_calc_filename = "xc_pot_calc.dat";

        string gphi_filename = "gphi.dat";

        string pot_result_filename = "y.dat";
        string laplacian_file = "lap.dat";

        double top_bc, split_bc1, split_bc2, bottom_bc;
        List<double> gate_bcs;
        double split_width;
        double surface;
        
        double z_scaling;

        public TwoD_PoissonSolver_Scaled(Experiment exp, bool using_external_code, Dictionary<string, object> input)
            : base(using_external_code, input)
        {
            this.exp = exp;
            // check if the top layer is air... if so, we need to use natural boundary conditions on the upper surface (which must be zero)
            natural_topbc = (exp.Layers[exp.Layers.Length - 1].Material == Material.Air);

            // calculate scaling factor w = a * z such that the z dimension has the same length as the y dimension
            z_scaling = (exp.Ny_Pot * exp.Dy_Pot) / (exp.Nz_Pot * exp.Dz_Pot);
        }

        protected override Band_Data Parse_Potential(string location, string[] data)
        {
            string[] new_data = Trim_Potential_File(data);
            Band_Data result = Band_Data.Parse_Band_Data(location, new_data, exp.Ny_Dens, exp.Nz_Dens);

            if (File.Exists(laplacian_file))
            {
                // also parse the laplacian data from file
                string[] lines = File.ReadAllLines(laplacian_file);
                string[] lap_data = Trim_Potential_File(lines);
                result.Laplacian = Band_Data.Parse_Band_Data(laplacian_file, lap_data, exp.Ny_Dens, exp.Nz_Dens);
            }
            else
                // if the file doesn't exist, initiate an empty laplacian
                result.Initiate_Laplacian();

            return result;
        }

        protected override void Save_Data(Band_Data density, string input_file_name)
        {
            // check that the data isn't all zeros and if it is, add a small 1e-8 perturbation at the centre
            if (density.mat.Min() > -1e-8 && density.mat.Max() < 1e-8)
                density.mat[(int)(density.mat.Rows / 2), (int)(density.mat.Cols / 2)] = 1e-8;

            density.Save_2D_Data(input_file_name, exp.Dy_Dens, z_scaling * exp.Dz_Dens, exp.Ymin_Dens, z_scaling * exp.Zmin_Dens);
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

            // and load generic gate boundary conditions
            gate_bcs = new List<double>();
            int count = 0;
            while (boundary_conditions.ContainsKey("V" + count.ToString()))
            { 
                gate_bcs.Add(boundary_conditions["V" + count.ToString()] * Physics_Base.energy_V_to_meVpzC);
                count++;
            }

            // and save the split width and surface charge
            this.split_width = device_dimension["split_width"];
            this.surface = boundary_conditions["surface"];

            if (flexpde_script != null)
                Create_FlexPDE_File(top_bc, split_bc1, split_bc2, split_width, surface, bottom_bc, flexpde_script);
        }

        protected override Band_Data Get_Pot_From_External(Band_Data density)
        {
            Save_Data(density, dens_filename);

            return Get_Data_From_External(initcalc_location, flexpde_options, initcalc_result_filename);
        }

        public void Create_FlexPDE_File(double top_bc, double split_bc1, double split_bc2, double split_width, double surface, double bottom_bc, string output_file)
        {
            StreamWriter sw = new StreamWriter(output_file);

            sw.WriteLine("TITLE \'Split Gate\'");
            sw.WriteLine("COORDINATES cartesian2");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            sw.WriteLine("\tERRLIM=" + pot_tol.ToString());
            sw.WriteLine("\tGRIDLIMIT=20");
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\tband_gap");
            sw.WriteLine();
            // and the tables for carrier and donor densities
            sw.WriteLine("\trho_carrier = TABLE(\'" + dens_filename + "\')");
            sw.WriteLine("\trho_dopent = TABLE(\'" + densdopent_filename + "\')");
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\tz_scaling = " + z_scaling.ToString());
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine("\tlz = " + (exp.Dz_Pot * exp.Nz_Pot).ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = " + surface.ToString() + " * z_scaling ");
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in meV zC^-1)");
            for (int i = 0; i < gate_bcs.Count; i++)
                sw.WriteLine("\tV" + i.ToString() + " = " + gate_bcs[i].ToString());
            sw.WriteLine("\ttop_V = " + top_bc.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString());
            sw.WriteLine("\tsplit_depth = 10 * z_scaling\t! depth of the split gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + (exp.Layers[1].Zmax - 5).ToString() + " * z_scaling");
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivity");
            sw.WriteLine("\teps");
            sw.WriteLine();
            // other physical parameters
            sw.WriteLine("\tq_e = " + Physics_Base.q_e.ToString() + "! charge of electron in zC");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: dx(eps * dx(u)) + z_scaling * dy(eps * z_scaling * dy(u)) = - (rho_carrier + rho_dopent) \t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");
            Draw_Domain(sw);

            // add a front commmand at the well depth to introduce a higher density of points at the interface
            sw.WriteLine("\tFRONT(y - well_depth, z_scaling * 20)");
            sw.WriteLine();

            // write out plotting routines
            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(- q_e * u + 0.5 * band_gap) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tCONTOUR(u) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tCONTOUR(rho_carrier + rho_dopent) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(u) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(rho_carrier) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(- q_e * u + 0.5 * band_gap) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tCONTOUR(u) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (0, " + (z_scaling * (exp.Zmin_Dens + (exp.Nz_Dens + 1) * exp.Dz_Dens)).ToString() + ") TO (0, " + (z_scaling * (exp.Zmin_Dens - 2.0 * exp.Dz_Dens)).ToString() + ")");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (-" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth) TO (" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth)");
            sw.WriteLine("\tCONTOUR(rho_carrier + rho_dopent) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tELEVATION(rho_carrier + rho_dopent) FROM (0, " + (z_scaling * (exp.Zmin_Dens + (exp.Nz_Dens + 1) * exp.Dz_Dens)).ToString() + ") TO (0, " + (z_scaling * (exp.Zmin_Dens - 2.0 * exp.Dz_Dens)).ToString() + ")");
            sw.WriteLine("\tELEVATION(-1.0e21*rho_carrier/q_e) FROM (-" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth) TO (" + Math.Abs(1.2 * exp.Ymin_Dens).ToString() + ", well_depth)");

            // and transfer the data to a file for reloading and replotting later
            sw.WriteLine();
            sw.WriteLine("\tTABLE(u) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + initcalc_result_filename + "\"");
            sw.WriteLine("\tTABLE(-1.0 * rho_carrier - rho_dopent) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + laplacian_file + "\"");
            sw.WriteLine("\tTRANSFER(u) FILE = \'" + pot_filename + "\'");
            sw.WriteLine("\tTRANSFER(0.0 * u) FILE = \'" + new_pot_filename + "\' ! dummy file for smoother function");
            sw.WriteLine();

            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }

        void Draw_Domain(StreamWriter sw)
        {
            // cycle through layers below surface
            int count = 1;
            int voltage_count = 0;
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                for (int j = 0; j < exp.Layers[i].No_Components; j++)
                {
                    ILayer current_layer = exp.Layers[i].Get_Component(j);

                    // set the minimum and maximum positions (which are limited here by the FlexPDE solution domain)
                    string ymin, ymax;
                    if (current_layer.Ymin != double.MinValue)
                        ymin = current_layer.Ymin.ToString();
                    else ymin = "-ly / 2";
                    if (current_layer.Ymax != double.MaxValue)
                        ymax = current_layer.Ymax.ToString();
                    else ymax = "ly / 2";

                    // output region data
                    sw.WriteLine("\tREGION " + count.ToString());
                    sw.WriteLine("\t\teps = " + current_layer.Permitivity.ToString());
                    sw.WriteLine("\t\tband_gap = " + current_layer.Band_Gap.ToString());

                    // start defining the domain
                    sw.WriteLine("\t\tSTART(" + ymax + ", " + exp.Layers[i].Zmin.ToString() + " * z_scaling)");
                    // if a metal, use a predefined voltage
                    if (current_layer.Material == Material.Metal)
                    {
                        sw.WriteLine("\t\tVALUE(u) = V" + voltage_count.ToString());
                        voltage_count++;
                    }
                    sw.WriteLine("\t\tLINE TO (" + ymax + ", " + exp.Layers[i].Zmax.ToString() + " * z_scaling)");

                    // set top gate here
                    if (i == exp.Layers.Length - 1)
                    {
                        // sw.WriteLine("\t\tVALUE(u) = split_V\n\t\tline TO (split_width / 2, 0)\n\t\tNATURAL(u) = surface_bc\n\t\tLINE TO (-split_width / 2, 0)\n\t\tVALUE(u) = split_V");
                        sw.WriteLine("\t\tVALUE(u) = top_V");
                        if (natural_topbc)
                            sw.WriteLine("\t\tNATURAL(u) = top_V");
                    }
                    // or surface condition
                    if (exp.Layers[i].Zmax == 0.0)
                        sw.WriteLine("\t\tNATURAL(u) = surface_bc * (ustep(x + split_width / 2 - 20) - ustep(x - split_width / 2 + 20))");

                    sw.WriteLine("\t\tLINE TO (" + ymin + ", " + exp.Layers[i].Zmax.ToString() + " * z_scaling)");
                    // reset the boundary coondition if this is the top or surface
                    if (i == exp.Layers.Length - 1 || exp.Layers[i].Zmax == 0.0)
                        sw.WriteLine("\t\tNATURAL(u) = 0");
                    sw.WriteLine("\t\tLINE TO (" + ymin + ", " + exp.Layers[i].Zmin.ToString() + " * z_scaling)");
                    // set bottom boundary condition
                    if (i == 1)
                        sw.WriteLine("\t\tVALUE(u) = bottom_bc");
                    sw.WriteLine("\t\tLINE TO CLOSE");
                    sw.WriteLine();

                    count++;
                }
            }

        }

        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi)
        {
            Save_Data(gphi, gphi_filename);
            Save_Data(rho_prime.Spin_Summed_Data, densderiv_filename);
            Create_NewtonStep_File(split_width, flexpde_script, T);

            return Get_Data_From_External(newton_location, flexpde_options, newton_result_filename);
        }
        
        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data car_dens, Band_Data dft_pot, Band_Data dft_calc)
        {
            Save_Data(car_dens.Spin_Summed_Data, dens_filename);
            Save_Data(dft_pot, xc_pot_filename);
            Save_Data(dft_calc, xc_pot_calc_filename);

            return Calculate_Newton_Step(rho_prime, gphi);
        }

        public void Create_NewtonStep_File(double split_width, string output_file, double t)
        {
            StreamWriter sw = new StreamWriter(output_file);

            // an edge to the domain which FlexPDE will see gphi outside of the grid on which the density, etc. is defined.
            // the point is that the step functions will cut off the edge of lap.dat on one side so we make the window a little bigger
            double gphi_tol = 1.0;

            // generate summary containing residual integrals for car dens (will be used later for x and phi)
            string window_function_string;
            if (exp.Ymin_Dens < 0.0)
                window_function_string = "(ustep(x + " + (-1.0 * (exp.Ymin_Dens - gphi_tol)).ToString() + ")";
            else
                window_function_string = "(ustep(x - " + (exp.Ymin_Dens - gphi_tol).ToString() + ")";
            double xmax = exp.Ymin_Dens + (exp.Ny_Dens - 1) * exp.Dy_Dens;
            if (xmax < 0.0)
                window_function_string = window_function_string + " - ustep(x + " + (-1.0 * (xmax + gphi_tol)).ToString() + "))";
            else
                window_function_string = window_function_string + " - ustep(x - " + (xmax + gphi_tol).ToString() + "))";
            if (exp.Zmin_Dens < 0.0)
                window_function_string = window_function_string + " * (ustep(y + " + (-1.0 * z_scaling * (exp.Zmin_Dens - gphi_tol)).ToString() + ")";
            else
                window_function_string = window_function_string + " * (ustep(y - " + (z_scaling * (exp.Zmin_Dens - gphi_tol)).ToString() + ")";
            double ymax = z_scaling * (exp.Zmin_Dens + (exp.Nz_Dens - 1) * exp.Dz_Dens);
            if (ymax < 0.0)
                window_function_string = window_function_string + " - ustep(y + " + (-1.0 * (ymax + gphi_tol)).ToString() + "))";
            else
                window_function_string = window_function_string + " - ustep(y - " + (ymax + gphi_tol).ToString() + "))";

            //string minus_g_phi = "(dx(eps * dx(phi + t * new_phi)) + z_scaling * dy(eps * z_scaling * dy(phi + t * new_phi)) + " +  window_function_string + " * car_dens + dop_dens)";
            //minus_g_phi += " * upulse(y - well_depth + 200, y - well_depth - 100)";
            string g_phi = "gphi * " + window_function_string;

            sw.WriteLine("TITLE \'Split Gate\'");
            sw.WriteLine("COORDINATES cartesian2");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            sw.WriteLine("\tERRLIM=" + newton_tol.ToString());
            // no regridding for the newton step.  Just use the original grid from the potential calculation
            //sw.WriteLine("REGRID = OFF");
            sw.WriteLine();
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\tband_gap");
            sw.WriteLine();
            sw.WriteLine("\tTRANSFERMESH(\'" + pot_filename + "\', phi_old)");
            sw.WriteLine("\tTRANSFER(\'" + new_pot_filename + "\', x_old)");
            sw.WriteLine();
            // and the tables for carrier and donor densities
            //sw.WriteLine("\tcar_dens = SMOOTH(" + exp.Dy_Dens.ToString() + ") TABLE(\'" + dens_filename + "\')");
            sw.WriteLine("\tgphi = SPLINE TABLE(\'" + gphi_filename + "\')");
            sw.WriteLine("\txc_pot = SPLINE TABLE(\'" + xc_pot_filename + "\')");
            sw.WriteLine("\txc_pot_calc = SPLINE TABLE(\'xc_pot_calc.dat\')");
            sw.WriteLine("\tcar_dens = SPLINE TABLE(\'" + dens_filename + "\')");
            sw.WriteLine("\tdop_dens = TABLE(\'" + densdopent_filename + "\')");
            sw.WriteLine("\trho_prime = TABLE(\'" + densderiv_filename + "\')");
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\tz_scaling = " + z_scaling.ToString());
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine("\tlz = " + (exp.Dz_Pot * exp.Nz_Pot).ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = 0.0");
            sw.WriteLine("\tsurface_bc = 0.0");
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in meV zC^-1)");
            for (int i = 0; i < gate_bcs.Count; i++)
                sw.WriteLine("\tV" + i.ToString() + " = 0.0");
            sw.WriteLine("\ttop_V = 0.0");
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString());
            sw.WriteLine("\tsplit_depth = 10 * z_scaling\t! scaled depth of the split gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + (exp.Layers[1].Zmax - 5).ToString() + " * z_scaling");
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
            sw.WriteLine("\tu: -1.0 * dx(eps * dx(u)) - 1.0 * z_scaling * dy(eps * z_scaling * dy(u)) - rho_prime * u = -1.0 * " + g_phi + " \t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the different layers
            sw.WriteLine("BOUNDARIES");
            Draw_Domain(sw);
            
            sw.WriteLine("\tRESOLVE(-1.0 * " + g_phi + ")");
            sw.WriteLine();

            // write out plotting routines

            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(-1.0 * " + g_phi + ") ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tCONTOUR(u) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tELEVATION(u) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(-1.0 * " + g_phi + ") FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(phi_old + t * x_old) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tCONTOUR((phi_old + t * x_old) * q_e) ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(car_dens * " + window_function_string + ") ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(-1.0 * " + g_phi + ") ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(u) ON GRID(x, y / z_scaling)");
            sw.WriteLine("\tCONTOUR(u * q_e) ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(xc_pot) ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(xc_pot_calc) ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            sw.WriteLine("\tCONTOUR(rho_prime) ON GRID(x, y / z_scaling) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ")");
            
            window_function_string = "(1.0 - " + window_function_string + ") - (ustep(y + 3) - ustep(y -13)) * (ustep(x + 353) - ustep(x + 347) - ustep(x - 353) + ustep(x - 347))";

//            string residual_report_filename = "residual_g.dat";
//            sw.WriteLine("\tSUMMARY EXPORT FILE = \'" + residual_report_filename + "\'");
//            sw.WriteLine("\t\tREPORT ( INTEGRAL (-1.0 * (" + minus_g_phi + ") * u * " + window_function_string + ") / z_scaling) AS \"residual_g_phi\"");
//            sw.WriteLine("\t\tREPORT ( INTEGRAL ((dx(eps * dx(u)) + z_scaling * dy(eps * z_scaling * dy(u))) * u * " + window_function_string + ") / z_scaling) AS \"residual_g_x\"");

            // and transfer the data to a file for reloading and replotting later
            sw.WriteLine();
            sw.WriteLine("\tTRANSFER(phi_old + t * x_old) FILE = \'" + pot_filename + "\'");
            sw.WriteLine("\tTRANSFER(u) FILE = \'" + new_pot_filename + "\'");
            sw.WriteLine("\tTABLE(phi_old + t * x_old) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + pot_result_filename + "\"");
            sw.WriteLine("\tTABLE(" + g_phi + " - rho_prime * u) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + laplacian_file + "\"");
            sw.WriteLine("\tTABLE(u) ZOOM (" + exp.Ymin_Dens.ToString() + ", " + (z_scaling * exp.Zmin_Dens).ToString() + ", " + ((exp.Ny_Dens - 1) * exp.Dy_Dens).ToString() + ", " + (z_scaling * (exp.Nz_Dens - 1) * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"" + newton_result_filename + "\"");
            sw.WriteLine();

            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }

        public override Band_Data Chemical_Potential
        {
            get
            {
                string[] lines = File.ReadAllLines(pot_result_filename);
                string[] data = Trim_Potential_File(lines);

                // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
                return Physics_Base.q_e * Parse_Potential(pot_result_filename, data);
            }
        }
    }
}
