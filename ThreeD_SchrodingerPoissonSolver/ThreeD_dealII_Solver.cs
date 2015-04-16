using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Geometry;
using CenterSpace.NMath.Core;

namespace ThreeD_SchrodingerPoissonSolver
{
    class ThreeD_dealII_Solver : dealII_Base
    {
        Experiment exp;

        public ThreeD_dealII_Solver(Experiment exp, bool using_external_code, Dictionary<string, object> input)
            : base(using_external_code, input)
        {
            this.exp = exp;
        }

        /// <summary>
        /// initiates the poisson solver by writing the input file for the parameter handlers of the initcalc and newton methods
        /// </summary>
        public override void Initiate_Poisson_Solver(Dictionary<string, double> device_dimensions, Dictionary<string, double> boundary_conditions)
        {
            // write the parameter details for the initial potential calculation
            StreamWriter sw_initcalc = new StreamWriter(initcalc_parameterfile);

            // check if the top layer is air... if so, we need to use natural boundary conditions on the upper surface (this is almost always going to be the case in 3D)
            bool natural_top_bc = (exp.Layers[exp.Layers.Length - 1].Material == Material.Air);
            if (natural_top_bc)
                device_dimensions["top_V"] = 0.0;

            sw_initcalc.WriteLine("subsection System Geometry");
            sw_initcalc.WriteLine("\tset zmin_pot = " + device_dimensions["zmin_pot"].ToString());
            sw_initcalc.WriteLine("\tset split_width = " + device_dimensions["split_width"].ToString());
            sw_initcalc.WriteLine("\tset split_length = " + device_dimensions["split_length"].ToString());
            sw_initcalc.WriteLine("\tset top_length = " + device_dimensions["top_length"].ToString());
            sw_initcalc.WriteLine("\tset pmma_depth = " + device_dimensions["pmma_depth"].ToString());
            sw_initcalc.WriteLine("\tset cap_depth = " + device_dimensions["cap_depth"].ToString());
            sw_initcalc.WriteLine("\tset interface_depth = " + device_dimensions["interface_depth"].ToString());
            sw_initcalc.WriteLine("\tset buffer_depth = " + device_dimensions["buffer_depth"].ToString());
            sw_initcalc.WriteLine("end");

            sw_initcalc.WriteLine("subsection Boundary Conditions");
            sw_initcalc.WriteLine("\tset top_bc = " + (boundary_conditions["top_V"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset bottom_bc = " + (boundary_conditions["bottom_V"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset split_bc1 = " + (boundary_conditions["split_V1"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset split_bc2 = " + (boundary_conditions["split_V2"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset surface_bc = " + (0.5 * boundary_conditions["surface"]).ToString());                                   // note that the surface "kink" is halved...
                                                                                                                                                // not sure why, but this is deal.II specific
            sw_initcalc.WriteLine("end");

            sw_initcalc.WriteLine("subsection Carrier Density Parameters");
            sw_initcalc.WriteLine("\tset nx_dens = " + exp.Nx_Dens.ToString());
            sw_initcalc.WriteLine("\tset ny_dens = " + exp.Ny_Dens.ToString());
            sw_initcalc.WriteLine("\tset nz_dens = " + exp.Nz_Dens.ToString());
            sw_initcalc.WriteLine("\tset xmin_dens = " + exp.Xmin_Dens.ToString());
            sw_initcalc.WriteLine("\tset xmax_dens = " + (exp.Xmin_Dens + exp.Dx_Dens * (exp.Nx_Dens - 1)).ToString());
            sw_initcalc.WriteLine("\tset ymin_dens = " + exp.Ymin_Dens.ToString());
            sw_initcalc.WriteLine("\tset ymax_dens = " + (exp.Ymin_Dens + exp.Dy_Dens * (exp.Ny_Dens - 1)).ToString());
            sw_initcalc.WriteLine("\tset zmin_dens = " + exp.Zmin_Dens.ToString());
            sw_initcalc.WriteLine("\tset zmax_dens = " + (exp.Zmin_Dens + exp.Dz_Dens * (exp.Nz_Dens - 1)).ToString());
            sw_initcalc.WriteLine("end");

            sw_initcalc.WriteLine("subsection Donor Density Parameters");
            sw_initcalc.WriteLine("\tset nz_donor = " + nz_donor);
            sw_initcalc.WriteLine("\tset zmin_donor = " + zmin_donor);
            sw_initcalc.WriteLine("\tset zmax_donor = " + zmax_donor);
            sw_initcalc.WriteLine("end");

            sw_initcalc.WriteLine("set natural top_bc = " + natural_top_bc.ToString().ToLower());

            sw_initcalc.Close();

            // and write the parameter details needed for the newton step
            StreamWriter sw_newton = new StreamWriter(newton_parameterfile);
            sw_newton.WriteLine("subsection System Geometry");
            sw_newton.WriteLine("\tset zmin_pot = " + device_dimensions["zmin_pot"].ToString());
            sw_newton.WriteLine("\tset split_width = " + device_dimensions["split_width"].ToString());
            sw_newton.WriteLine("\tset split_length = " + device_dimensions["split_length"].ToString());
            sw_newton.WriteLine("\tset top_length = " + device_dimensions["top_length"].ToString());
            sw_newton.WriteLine("\tset pmma_depth = " + device_dimensions["pmma_depth"].ToString());
            sw_newton.WriteLine("\tset cap_depth = " + device_dimensions["cap_depth"].ToString());
            sw_newton.WriteLine("\tset interface_depth = " + device_dimensions["interface_depth"].ToString());
            sw_newton.WriteLine("\tset buffer_depth = " + device_dimensions["buffer_depth"].ToString());
            sw_newton.WriteLine("end");

            sw_newton.WriteLine("subsection Carrier Density Parameters");
            sw_newton.WriteLine("\tset nx_dens = " + exp.Nx_Dens.ToString());
            sw_newton.WriteLine("\tset ny_dens = " + exp.Ny_Dens.ToString());
            sw_newton.WriteLine("\tset nz_dens = " + exp.Nz_Dens.ToString());
            sw_newton.WriteLine("\tset xmin_dens = " + exp.Xmin_Dens.ToString());
            sw_newton.WriteLine("\tset xmax_dens = " + (exp.Xmin_Dens + exp.Dx_Dens * (exp.Nx_Dens - 1)).ToString());
            sw_newton.WriteLine("\tset ymin_dens = " + exp.Ymin_Dens.ToString());
            sw_newton.WriteLine("\tset ymax_dens = " + (exp.Ymin_Dens + exp.Dy_Dens * (exp.Ny_Dens - 1)).ToString());
            sw_newton.WriteLine("\tset zmin_dens = " + exp.Zmin_Dens.ToString());
            sw_newton.WriteLine("\tset zmax_dens = " + (exp.Zmin_Dens + exp.Dz_Dens * (exp.Nz_Dens - 1)).ToString());
            sw_newton.WriteLine("end");

            sw_newton.WriteLine("set natural top_bc = " + natural_top_bc.ToString().ToLower());

            sw_newton.Close();
        }

        protected override Band_Data Parse_Potential(string[] data)
        {
            string[] new_data = Trim_Potential_File(data);
            return Band_Data.Parse_Band_Data(new_data, exp.Nx_Dens, exp.Ny_Dens, exp.Nz_Dens);
        }

        protected override void Save_Data(Band_Data density, string input_file_name)
        {
            density.Save_Data(input_file_name);
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
    }
}
