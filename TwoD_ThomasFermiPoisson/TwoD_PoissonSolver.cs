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
    public class TwoD_PoissonSolver : Potential_Base
    {
        double bottom_bc;
        Experiment exp;

        public TwoD_PoissonSolver(Experiment exp, bool using_flexPDE, string flexPDE_input, string flexPDE_location, double tol)
            : base(using_flexPDE, flexPDE_input, flexPDE_location, tol)
        {
            this.exp = exp;
            this.dens_filename = "dens_2D.dat";
        }

        protected override Band_Data Parse_Potential(string[] data)
        {
            // and check that there is the right number of data points back
            if (data.Length != exp.Ny_Dens * exp.Nz_Dens)
                throw new Exception("Error - FlexPDE is outputting the wrong number of potential data points");

            // and parse these values into a DoubleVector
            Band_Data result = new Band_Data(new DoubleMatrix(exp.Ny_Dens, exp.Nz_Dens));
            for (int i = 0; i < exp.Ny_Dens; i++)
            {
                for (int j = 0; j < exp.Nz_Dens; j++)
                    result.mat[i, j] = double.Parse(data[j * exp.Ny_Dens + i]);
            }

            return result;
        }

        public void Create_FlexPDE_File(double surface, double bottom_bc, string output_file)
        {
            StreamWriter sw = new StreamWriter(output_file);

            sw.WriteLine("TITLE \'Split Gate\'");
            sw.WriteLine("COORDINATES cartesian2");
            sw.WriteLine("VARIABLES");
            sw.WriteLine("\tu");
            sw.WriteLine("SELECT");
            // gives the flexPDE tolerance for the finite element solve
            sw.WriteLine("\tERRLIM=1e-5");
            sw.WriteLine("DEFINITIONS");
            // this is where the density variable
            sw.WriteLine("\trho");
            sw.WriteLine("\tband_gap");
            sw.WriteLine();
            // simulation dimension
            sw.WriteLine("\tly = " + (exp.Dy_Pot * exp.Ny_Pot).ToString());
            sw.WriteLine("\tlz = " + (exp.Dz_Pot * exp.Nz_Pot).ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = " + surface.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in V)");
            sw.WriteLine("\tsplit_V = " + split_V.ToString());
            sw.WriteLine("\ttop_V = " + top_V.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = 600");
            sw.WriteLine("\tsplit_depth = 10\t! depth of the split gate metal material");
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
            // other physical parameters
            sw.WriteLine("\tq_e = " + Physics_Base.q_e.ToString() + "! charge of electron in zC");
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: div(eps * grad(u)) = -rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");

            // cycle through layers below surface
            for (int i = 1; i < exp.Layers.Length; i++)
            {
                sw.WriteLine("\tREGION " + exp.Layers[i].Layer_No.ToString());
                if (exp.Layers[i].Acceptor_Conc != 0.0 || exp.Layers[i].Donor_Conc != 0.0 || exp.Layers[i].Layer_No == Geom_Tool.Find_Layer_Below_Surface(exp.Layers))
                    sw.WriteLine("\t\trho = TABLE(\'dens_2D_donors.dat\', x, y)");
                else if (exp.Layers[i].Layer_No <= Geom_Tool.Find_Layer_Below_Surface(exp.Layers))
                    sw.WriteLine("\t\trho = TABLE(\'dens_2D.dat\', x, y)");
                else
                    sw.WriteLine("\t\trho = 0.0");
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(exp.Layers[i].Material));
                sw.WriteLine("\t\tband_gap = " + exp.Layers[i].Band_Gap.ToString());
                sw.WriteLine("\t\tSTART(ly / 2, " + exp.Layers[i].Zmin.ToString() + ")");
                sw.WriteLine("\t\tLINE TO (ly / 2, " + exp.Layers[i].Zmax.ToString() + ")");
                // set top gate here
                if (i == exp.Layers.Length - 1)
                    sw.WriteLine("\t\tVALUE(u) = top_V");
                // or surface condition
                if (i == Geom_Tool.Find_Layer_Below_Surface(exp.Layers))
                {
                    sw.WriteLine("\t\tLINE TO (-split_width / 2, " + exp.Layers[i].Zmax.ToString() + ")");
                    sw.WriteLine("\t\tNATURAL(u) = surface_bc");
                    sw.WriteLine("\t\tLINE TO (split_width / 2, " + exp.Layers[i].Zmax.ToString() + ")");
                }
                sw.WriteLine("\t\tLINE TO (-ly / 2, " + exp.Layers[i].Zmax.ToString() + ")");
                sw.WriteLine("\t\tNATURAL(u) = 0 LINE TO (-ly / 2, " + exp.Layers[i].Zmin.ToString() + ")");
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

            // write out plotting routines

            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tCONTOUR(-q_e * u + 0.5 * band_gap)");

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tELEVATION(-q_e * u + 0.5 * band_gap) FROM (0, 0) TO (0, -lz)");
            sw.WriteLine("\tELEVATION(-q_e * u + 0.5 * band_gap) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");
            sw.WriteLine("\tCONTOUR(-q_e * u + 0.5 * band_gap)");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tELEVATION(rho) FROM (0, 0) TO (0, -lz)");
            sw.WriteLine("\tELEVATION(rho) FROM (-ly / 2, well_depth) TO (ly / 2, well_depth)");

            // and transfer the data to a file for reloading and replotting later
            sw.WriteLine();
            sw.WriteLine("\tTABLE(u) ZOOM ("+ exp.Ymin_Dens.ToString() + ", " + exp.Zmin_Dens.ToString() + ", " + (exp.Ny_Dens * exp.Dy_Dens).ToString() + ", " + (exp.Nz_Dens * exp.Dz_Dens).ToString() + ") EXPORT FORMAT \"#1\" POINTS = (" + exp.Ny_Dens.ToString() + ", " + exp.Nz_Dens.ToString() + ") FILE = \"pot.dat\"");
            sw.WriteLine("\tTRANSFER (rho, u, -q_e * u + 0.5 * band_gap) FILE=\"data_file.dat\"");
            sw.WriteLine();

            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }

        double top_V, split_V;
        public void Set_Boundary_Conditions(double top_V, double split_V, double bottom_bc, double surface)
        {
            // change the boundary conditions to potential boundary conditions by dividing through by -q_e
            // (as phi = E_c / (-1.0 * q_e)
            this.top_V = top_V; this.split_V = split_V;
            this.bottom_bc = bottom_bc / (-1.0 * Physics_Base.q_e);

            if (flexpde_inputfile != null)
                Create_FlexPDE_File(surface, bottom_bc, flexpde_inputfile);
        }

        protected override Band_Data Get_BandEnergy_On_Regular_Grid(Band_Data density)
        {
            throw new NotImplementedException();
        }

        protected override void Save_Density_Data(Band_Data density, string input_file_name)
        {
            density.Save_2D_Data(input_file_name, exp.Dy_Dens, exp.Dz_Dens, exp.Ymin_Dens, exp.Zmin_Dens);
        }
    }
}
