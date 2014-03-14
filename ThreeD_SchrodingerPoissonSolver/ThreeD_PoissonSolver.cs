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
    public class ThreeD_PoissonSolver: Potential_Base
    {
        double bottom_bc;
        Experiment exp;
        DoubleVector dens_1d;

        public ThreeD_PoissonSolver(Experiment exp, DoubleVector dens_1d, bool using_flexPDE, string flexPDE_input, string flexPDE_location, double tol)
            : base(using_flexPDE, flexPDE_input, flexPDE_location, tol)
        {
            this.exp = exp;
            this.dens_1d = dens_1d;
            this.dens_filename = "dens_3D.dat";
        }

        protected override Band_Data Parse_Potential(string[] data)
        {
            // and check that there is the right number of data points back
            if (data.Length != exp.Nx_Dens * exp.Ny_Dens)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleMatrix(exp.Nx_Dens, exp.Ny_Dens));
            for (int i = 0; i < exp.Nx_Dens; i++)
            {
                for (int j = 0; j < exp.Ny_Dens; j++)
                    result.mat[i, j] = double.Parse(data[j * exp.Nx_Dens + i]);
            }

            return result;
        }

        public void Create_FlexPDE_File(double split_width, double split_length, double surface, double bottom_bc, string output_file)
        {
            StreamWriter sw = new StreamWriter(output_file);

            // write out output file
            sw.WriteLine("TITLE \'Full Split Gate Geometry\'");
            sw.WriteLine("COORDINATES cartesian3");
            sw.WriteLine();
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine();
            sw.WriteLine("SELECT");
            sw.WriteLine("\tERRLIM=1e-3");
            sw.WriteLine("DEFINITIONS");
            sw.WriteLine("\trho = 0.0");
            sw.WriteLine("\tband_gap");
            sw.WriteLine();
            sw.WriteLine("\tlx = " + (exp.Dx_Pot * exp.Nx_Pot).ToString());
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine();
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = " + surface.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in V)");
            sw.WriteLine("\tsplit_V = " + split_V.ToString());
            sw.WriteLine("\ttop_V = " + top_V.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = " + split_width.ToString());
            sw.WriteLine("\tsplit_length = " + split_length.ToString());
            sw.WriteLine("\tsplit_depth = 10 ! depth of the split gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = " + (exp.Layers[1].Zmax - 1).ToString());
            sw.WriteLine();
            sw.WriteLine("\t! Electrical permitivities");
            sw.WriteLine("\teps_0 = " + Physics_Base.epsilon_0.ToString());
            // relative permitivity of materials
            sw.WriteLine("\teps_r_GaAs = " + Physics_Base.epsilon_r_GaAs.ToString());
            sw.WriteLine("\teps_r_AlGaAs = " + Physics_Base.epsilon_r_AlGaAs.ToString());
            sw.WriteLine("\teps_pmma = " + Physics_Base.epsilon_pmma.ToString());
            sw.WriteLine("\teps");
            sw.WriteLine();
            sw.WriteLine("\tq_e = " + Physics_Base.q_e.ToString() + "! charge of electron in zC");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            sw.WriteLine("\tu: div(eps * grad(u)) = - rho	! Poisson's equation");
            sw.WriteLine();
            sw.WriteLine("EXTRUSION");
            sw.WriteLine("\tSURFACE \"Substrate\"\tz = " + exp.Layers[0].Zmax.ToString());
            int layercount = 1;
            for (int i = 1; i < exp.Layers.Length + 1; i++)
            {
                sw.WriteLine("\t\tLAYER \"" + i.ToString() + "\"");
                if (i == Geom_Tool.Find_Layer_Above_Surface(exp.Layers))
                {
                    sw.WriteLine("\tSURFACE \"" + i.ToString() + "\"\tz = split_depth");
                    layercount--;
                }
                else
                    sw.WriteLine("\tSURFACE	\"" + i.ToString() + "\"\tz = " + exp.Layers[layercount].Zmax.ToString());

                layercount++;
            }
            sw.WriteLine();
            sw.WriteLine("BOUNDARIES");
            sw.WriteLine("\tSURFACE \"Substrate\"	VALUE(u) = bottom_bc");
            sw.WriteLine("\tSURFACE \"" + Geom_Tool.Find_Layer_Below_Surface(exp.Layers).ToString() + "\" NATURAL(u) = surface_bc");
            sw.WriteLine("\tSURFACE \"" + exp.Layers.Length.ToString() + "\" VALUE(u) = top_V");
            sw.WriteLine();
            sw.WriteLine("\tREGION 1");
            layercount = 1;
            for (int i = 1; i < exp.Layers.Length + 1; i++)
            {
                sw.WriteLine("\t\tLAYER \"" + i.ToString() + "\"");
                if (exp.Layers[layercount].Acceptor_Conc != 0.0 || exp.Layers[layercount].Donor_Conc != 0.0 || exp.Layers[layercount].Layer_No == Geom_Tool.Find_Layer_Below_Surface(exp.Layers))
                    sw.WriteLine("\t\trho = TABLE(\'dens_3D_donors.dat\', x, y)");
                else if (exp.Layers[layercount].Layer_No <= Geom_Tool.Find_Layer_Below_Surface(exp.Layers))
                    sw.WriteLine("\t\trho = TABLE(\'dens_3D.dat\', x, y)");
                else
                    sw.WriteLine("\t\trho = 0.0");
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(exp.Layers[layercount].Material));
                sw.WriteLine("\t\tband_gap = " + exp.Layers[layercount].Band_Gap.ToString());
                sw.WriteLine();

                if (i == Geom_Tool.Find_Layer_Above_Surface(exp.Layers))
                    layercount--;

                layercount++;
            }
            sw.WriteLine();
            sw.WriteLine("\t\tSTART(-lx / 2, -ly / 2)");
            sw.WriteLine("\t\tLINE TO (-lx / 2, ly / 2)");
            sw.WriteLine("\t\tLINE TO (lx / 2, ly / 2)");
            sw.WriteLine("\t\tLINE TO (lx / 2, -ly / 2)");
            sw.WriteLine("\t\tLINE TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 2");
            sw.WriteLine("\t\tSURFACE \"" + Geom_Tool.Find_Layer_Below_Surface(exp.Layers).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tSURFACE \"" + Geom_Tool.Find_Layer_Above_Surface(exp.Layers).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tLAYER \"" + Geom_Tool.Find_Layer_Above_Surface(exp.Layers).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (-lx / 2, -split_length / 2)");
            sw.WriteLine("\t\tLINE TO (-lx / 2, split_length / 2)");
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (-split_width / 2, split_length / 2) TO (-split_width / 2, -split_length / 2) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tLIMITED REGION 3");
            sw.WriteLine("\t\tSURFACE \"" + Geom_Tool.Find_Layer_Below_Surface(exp.Layers).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tSURFACE \"" + Geom_Tool.Find_Layer_Above_Surface(exp.Layers).ToString() + "\" VALUE(u) = split_V");
            sw.WriteLine("\t\tLAYER \"" + Geom_Tool.Find_Layer_Above_Surface(exp.Layers).ToString() + "\" VOID");
            sw.WriteLine("\t\tSTART (lx / 2, -split_length / 2)");
            sw.WriteLine("\t\tLINE TO (lx / 2, split_length / 2)");
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (split_width / 2, split_length / 2) TO (split_width / 2, -split_length / 2) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(rho) ON z = well_depth");
            sw.WriteLine("\tCONTOUR(u) ON z = well_depth");
            sw.WriteLine("\tCONTOUR(u) ON x = 0");
            sw.WriteLine("\tCONTOUR(u) ON y = 0");
            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(rho) ON x = 0");
            sw.WriteLine("\tCONTOUR(u) ON x = 0");
            sw.WriteLine("\tCONTOUR(u) ON y = 0");
            sw.WriteLine("\tCONTOUR(rho) ON z = well_depth");
            sw.WriteLine("\tCONTOUR(- q_e * u + 0.5 * band_gap) ON z = well_depth");
            sw.WriteLine("\tELEVATION(rho) FROM (0,0, " + Geom_Tool.Get_Zmin(exp.Layers).ToString() + ") TO (0, 0, " + exp.Layers[exp.Layers.Length - 1].Zmax.ToString() + ")");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (0, 0, " + Geom_Tool.Get_Zmin(exp.Layers).ToString() + ") TO (0, 0, " + exp.Layers[exp.Layers.Length - 1].Zmax.ToString() + ")");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (0, -ly / 2, well_depth) TO (0, ly / 2, well_depth)");
            sw.WriteLine("\tELEVATION(- q_e * u + 0.5 * band_gap) FROM (-lx/2, 0, well_depth) TO (lx / 2, 0, well_depth)");
            sw.WriteLine();
            sw.WriteLine("\tCONTOUR(u) ON z = well_depth EXPORT FORMAT \"#1\" POINTS = (" + exp.Nx_Dens.ToString() + ", " + exp.Ny_Dens.ToString() + ") FILE = \"pot.dat\"");
            //sw.WriteLine("TRANSFER (rho, u) FILE=\"data_file.dat\"");
            sw.WriteLine();
            sw.WriteLine("END");

            sw.Close();
        }

        double top_V, split_V;
        public void Set_Boundary_Conditions(double top_V, double split_V, double split_width, double split_length, double bottom_bc, double surface)
        {
            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // (as phi = E_c / (-1.0 * q_e)
            this.top_V = top_V; this.split_V = split_V;
            this.bottom_bc = bottom_bc / (-1.0 * Physics_Base.q_e);

            if (flexpde_inputfile != null)
                Create_FlexPDE_File(split_width, split_length, surface, bottom_bc, flexpde_inputfile);
        }

        protected override Band_Data Get_ChemPot_On_Regular_Grid(Band_Data density)
        {
            throw new NotImplementedException();
        }

        protected override void Save_Density_Data(Band_Data density, string input_file_name)
        {
            density.Save_3D_Data(input_file_name, dens_1d, exp.Dx_Dens, exp.Dy_Dens, exp.Dz_Dens, exp.Xmin_Dens, exp.Ymin_Dens, exp.Zmin_Dens);
        }
    }
}
