using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace Split_Gate_Potential_Calc
{
    class SG_File_Generator
    {
        Dictionary<string, object> input;
        ILayer[] layers;
        double ly, lz;
        double bottom_bc;
        string rbf_fit;

        public SG_File_Generator(Dictionary<string, object> input)
        {
            this.input = input;
            this.layers = (ILayer[])input["Layers"];

            this.ly = (double)input["dy"] * (double)input["ny"];
            this.lz = (double)input["dz"] * (double)input["nz"];

            // Calculate the bottom boundary condition
            OneD_ThomasFermiPoisson.OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiPoisson.OneD_ThomasFermiSolver((double)input["T"], (double)input["dz"], (int)(double)input["nz"], (double)input["zmin"]);
            bottom_bc = dens_solv.Get_Chemical_Potential(Geom_Tool.Get_Zmin(layers), layers);

            // get the donor densities and positions
            Fit_Donor_Level_RBF(input);
        }

        void Fit_Donor_Level_RBF(Dictionary<string, object> input)
        {
            // find the donor layer
            int[] donor_layer = Layer_Tool.Get_Donor_Layers(layers);
            if (donor_layer.Length != 1) throw new NotImplementedException("Error - at the moment there can only be one donor layer in the sample!");

            // and get some of its parameters
            double donor_layer_width = layers[donor_layer[0]].Zmax - layers[donor_layer[0]].Zmin;
            int no_pos = (int)Math.Round(donor_layer_width / (double)input["dz"]);
            int donor_base_index = (int)(Math.Round(layers[donor_layer[0]].Zmin - Geom_Tool.Get_Zmin(layers)) / (double)input["dz"]);

            // calculate the positions and densities to make the RBF fit to
            double[] positions = new double[no_pos];
            double[] donor_dens = new double[no_pos];
            for (int i = 0; i < no_pos; i++)
            {
                positions[i] = layers[donor_layer[0]].Zmin + i * (double)input["dz"];
                donor_dens[i] = ((CenterSpace.NMath.Core.DoubleVector)input["oned_dens_data"])[donor_base_index + i];
            }

            // and do the fit here
            OneD_ThomasFermiPoisson.OneD_RBF_Fit rbf_fitter = new OneD_ThomasFermiPoisson.OneD_RBF_Fit(positions, 3.0 * (double)input["dz"]);
            rbf_fitter.Get_RBF_Weights(donor_dens);
            rbf_fit = rbf_fitter.Get_RBF_Equation("donor_dens", "y");
        }

        internal void Generate_2D_FlexPDE_File(double surface)
        {
            string output_file = "split_gate_2d.pde";

            if (File.Exists(output_file))
            {
                bool problem_resolved = false;
                while (!problem_resolved)
                {
                    Console.WriteLine(output_file + " already exists, do you want to replace it y/n");
                    char input_key = Console.ReadKey().KeyChar;
                    if (input_key == 'Y' || input_key == 'y')
                    {
                        File.Delete(output_file);
                        problem_resolved = true;
                    }
                    else if (input_key == 'N' || input_key == 'n')
                    {
                        Console.WriteLine("Input file name to save FlexPDE file to");
                        output_file = Console.ReadLine();
                        problem_resolved = true;
                    }
                    else
                        Console.WriteLine(input_key + " - Invalid input");
                }
            }

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
            sw.WriteLine("\tly = " + ly.ToString());
            sw.WriteLine("\tlz = " + lz.ToString());
            sw.WriteLine();
            // boundary conditions
            sw.WriteLine("\tbottom_bc = " + bottom_bc.ToString());
            sw.WriteLine("\tsurface_bc = " + surface.ToString());
            sw.WriteLine();
            sw.WriteLine("\t! GATE VOLTAGE INPUTS (in V)");
            sw.WriteLine("\tsplit_V = 0");
            sw.WriteLine("\ttop_V = 0");
            sw.WriteLine();
            sw.WriteLine("\t! SPLIT GATE DIMENSIONS (in nm)");
            sw.WriteLine("\tsplit_width = 600");
            sw.WriteLine("\tsplit_depth = 10\t! depth of the split gate metal material");
            sw.WriteLine();
            sw.WriteLine("\t! WELL DEPTH (in nm)");
            sw.WriteLine("\twell_depth = 100");
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
            // and the donor density functional fit
            sw.WriteLine(rbf_fit);
            sw.WriteLine();
            sw.WriteLine("EQUATIONS");
            // Poisson's equation
            sw.WriteLine("\tu: div(eps * grad(u)) = -rho\t! Poisson's equation");
            sw.WriteLine();
            // the boundary definitions for the differnet layers
            sw.WriteLine("BOUNDARIES");

            // cycle through layers below surface
            for (int i = 1; i < layers.Length; i++)
            {
                sw.WriteLine("\tREGION " + layers[i].Layer_No.ToString());
                if (layers[i].Acceptor_Conc != 0.0 || layers[i].Donor_Conc != 0.0)
                    //sw.WriteLine("\t\trho = TABLE(\'dens_2D.dat\', x, y)");
                    sw.WriteLine("\t\trho = donor_dens");
                else
                    sw.WriteLine("\t\trho = 0.0");
                sw.WriteLine("\t\teps = " + Layer_Tool.Get_Permitivity(layers[i].Material));
                sw.WriteLine("\t\tband_gap = " + layers[i].Band_Gap.ToString());
                sw.WriteLine("\t\tSTART(ly / 2, " + layers[i].Zmin.ToString() + ")");
                sw.WriteLine("\t\tLINE TO (ly / 2, " + layers[i].Zmax.ToString() + ")");
                // set top gate here
                if (i == layers.Length - 1)
                    sw.WriteLine("\t\tVALUE(u) = top_V");
                // or surface condition
                if (i == Geom_Tool.Find_Layer_Below_Surface(layers))
                {
                    sw.WriteLine("\t\tLINE TO (-split_width / 2, " + layers[i].Zmax.ToString() + ")");
                    sw.WriteLine("\t\tNATURAL(u) = surface_bc");
                    sw.WriteLine("\t\tLINE TO (split_width / 2, " + layers[i].Zmax.ToString() + ")");
                }
                sw.WriteLine("\t\tLINE TO (-ly / 2, " + layers[i].Zmax.ToString() + ")");
                sw.WriteLine("\t\tNATURAL(u) = 0 LINE TO (-ly / 2, " + layers[i].Zmin.ToString() + ")");
                // set bottom boundary condition
                if (i == 1)
                    sw.WriteLine("\t\tVALUE(u) = bottom_bc");
                sw.WriteLine("\t\tLINE TO CLOSE");
                sw.WriteLine();
            }

            // write in surface and gates
            sw.WriteLine("\tREGION " + layers.Length.ToString() + " ! Left split gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\tband_gap = 0");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(-ly / 2, 0)");
            // left split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (-ly / 2, split_depth) TO (-split_width / 2, split_depth) TO (-split_width / 2, 0) TO CLOSE");
            sw.WriteLine();
            sw.WriteLine("\tREGION " + (layers.Length + 1).ToString() + "! Right split gate");
            sw.WriteLine("\t\trho = 0");
            sw.WriteLine("\t\tband_gap = 0");
            sw.WriteLine("\t\teps = eps_0");
            sw.WriteLine("\t\tSTART(ly / 2, 0)");
            // right split gate voltage
            sw.WriteLine("\t\tVALUE(u) = split_V");
            sw.WriteLine("\t\tLINE TO (ly / 2, split_depth) TO (split_width / 2, split_depth) TO (split_width / 2, 0) TO CLOSE");
            sw.WriteLine();
            //sw.WriteLine("\tREGION " + (layers.Length + 2).ToString() + "! Surface charges");
            //sw.WriteLine("\t\trho = TABLE(\'dens_surface.dat\', x, y)");
            //sw.WriteLine("\t\tband_gap = " + layers[Geom_Tool.Find_Layer_Below_Surface(layers)].Band_Gap.ToString());
            //sw.WriteLine("\t\teps = eps_0");
            //sw.WriteLine("\t\tSTART(-split_width / 2, -0.1) LINE TO (-split_width / 2, 0.1) TO (split_width / 2, 0.1) TO (split_width / 2, -0.1) TO CLOSE");
            //sw.WriteLine();

            // write out plotting routines

            sw.WriteLine("MONITORS");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tCONTOUR(-q_e * u + 0.5 * band_gap)");

            sw.WriteLine("PLOTS");
            sw.WriteLine("\tCONTOUR(rho)");
            sw.WriteLine("\tCONTOUR(-q_e * u + 0.5 * band_gap)");
            sw.WriteLine("\tELEVATION(-q_e * u + 0.5 * band_gap) FROM (-ly / 2, -well_depth) TO (ly / 2, -well_depth)");
            sw.WriteLine("\tELEVATION(-q_e * u + 0.5 * band_gap) FROM (0, 0) TO (0, -lz)");
            sw.WriteLine("\tELEVATION(rho) FROM (0, 0) TO (0, -lz)");

            // and transfer the data to a file for reloading and replotting later
            sw.WriteLine();
            sw.WriteLine("\tTRANSFER (rho, u, -q_e * u + 0.5 * band_gap) FILE=\"data_file.dat\"");
            sw.WriteLine();

            sw.WriteLine("END");

            // and close the file writer
            sw.Close();
        }
    }
}
