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
using Solver_Bases;
using Solver_Bases.Geometry;
using CenterSpace.NMath.Core;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_dealII_Solver : dealII_Base
    {
        Experiment exp;
        string laplacian_file = "lap.dat";

        public TwoD_dealII_Solver(Experiment exp, bool using_external_code, Dictionary<string, object> input)
            : base(using_external_code, input)
        {
            this.exp = exp;
        }

        /// <summary>
        /// initiates the poisson solver by writing the input file for the parameter handlers of the initcalc and newton methods
        /// </summary>
        public override void Initiate_Poisson_Solver(Dictionary<string, double> device_dimensions, Dictionary<string, double> boundary_conditions)
        {
            // get split gate specific device dimensions
            device_dimensions["zmin_pot"] = exp.Layers[1].Zmin;
            device_dimensions["pmma_depth"] = exp.Layers[exp.Layers.Length - 1].Zmax;
            device_dimensions["cap_depth"] = Geom_Tool.Find_Layer_Below_Surface(exp.Layers).Zmin;
            device_dimensions["interface_depth"] = exp.Layers[1].Zmax;
            device_dimensions["buffer_depth"] = exp.Layers[2].Zmax;

            // write the parameter details for the initial potential calculation
            StreamWriter sw_initcalc = new StreamWriter(initcalc_parameterfile);

            // check if the top layer is air... if so, we need to use natural boundary conditions on the upper surface
            bool natural_top_bc = (exp.Layers[exp.Layers.Length - 1].Material == Material.Air);
            if (natural_top_bc)
                device_dimensions["top_V"] = 0.0;

            sw_initcalc.WriteLine("subsection System Geometry");
            sw_initcalc.WriteLine("\tset zmin_pot = " + device_dimensions["zmin_pot"].ToString());
            sw_initcalc.WriteLine("\tset split_width = " + device_dimensions["split_width"].ToString());
            sw_initcalc.WriteLine("\tset pmma_depth = " + device_dimensions["pmma_depth"].ToString());
            sw_initcalc.WriteLine("\tset cap_depth = " + device_dimensions["cap_depth"].ToString());
            sw_initcalc.WriteLine("\tset interface_depth = " + device_dimensions["interface_depth"].ToString());
            sw_initcalc.WriteLine("\tset buffer_depth = " + device_dimensions["buffer_depth"].ToString());
            sw_initcalc.WriteLine("end");

            sw_initcalc.WriteLine("subsection Boundary Conditions");
            sw_initcalc.WriteLine("\tset top_bc = " + (boundary_conditions["top_V"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset bottom_bc = " + (boundary_conditions["bottom_V"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset split_bc1 = " + (boundary_conditions["V0"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset split_bc2 = " + (boundary_conditions["V1"] * Physics_Base.energy_V_to_meVpzC).ToString());
            sw_initcalc.WriteLine("\tset surface_bc = " + (0.5 * boundary_conditions["surface"]).ToString());                                   // note that the surface "kink" is halved...
            // not sure why, but this is deal.II specific
            sw_initcalc.WriteLine("end");

            sw_initcalc.WriteLine("subsection Carrier Density Parameters");
            sw_initcalc.WriteLine("\tset ny_dens = " + exp.Ny_Dens.ToString());
            sw_initcalc.WriteLine("\tset nz_dens = " + exp.Nz_Dens.ToString());
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
            sw_newton.WriteLine("\tset pmma_depth = " + device_dimensions["pmma_depth"].ToString());
            sw_newton.WriteLine("\tset cap_depth = " + device_dimensions["cap_depth"].ToString());
            sw_newton.WriteLine("\tset interface_depth = " + device_dimensions["interface_depth"].ToString());
            sw_newton.WriteLine("\tset buffer_depth = " + device_dimensions["buffer_depth"].ToString());
            sw_newton.WriteLine("end");

            sw_newton.WriteLine("subsection Carrier Density Parameters");
            sw_newton.WriteLine("\tset ny_dens = " + exp.Ny_Dens.ToString());
            sw_newton.WriteLine("\tset nz_dens = " + exp.Nz_Dens.ToString());
            sw_newton.WriteLine("\tset ymin_dens = " + exp.Ymin_Dens.ToString());
            sw_newton.WriteLine("\tset ymax_dens = " + (exp.Ymin_Dens + exp.Dy_Dens * (exp.Ny_Dens - 1)).ToString());
            sw_newton.WriteLine("\tset zmin_dens = " + exp.Zmin_Dens.ToString());
            sw_newton.WriteLine("\tset zmax_dens = " + (exp.Zmin_Dens + exp.Dz_Dens * (exp.Nz_Dens - 1)).ToString());
            sw_newton.WriteLine("end");

            sw_newton.WriteLine("set natural top_bc = " + natural_top_bc.ToString().ToLower());

            sw_newton.Close();
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
    }
}
