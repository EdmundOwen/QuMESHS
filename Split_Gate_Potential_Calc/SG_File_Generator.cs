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
        int ny, nz;
        double dy, dz;
        double ymin, zmin;
        double bottom_bc;
        string rbf_fit;

        public SG_File_Generator(Dictionary<string, object> input)
        {
            this.input = input;
            this.layers = (ILayer[])input["Layers"];

            this.dy = (double)input["dy"]; this.ny = (int)(double)input["ny"];
            this.dz = (double)input["dz"]; this.nz = (int)(double)input["nz"];
            this.ymin = (double)input["ymin"]; this.zmin = (double)input["zmin"];

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

        internal void Generate_2D_FlexPDE_File(TwoD_ThomasFermiPoisson.Experiment exp, double surface)
        {
            string output_file = "split_gate_2d.pde";
            
            Check_Output_Filename(ref output_file);

            TwoD_ThomasFermiPoisson.TwoD_PoissonSolver temp_poissolv = new TwoD_ThomasFermiPoisson.TwoD_PoissonSolver(exp, true, output_file, null, 0.0);
            temp_poissolv.Create_FlexPDE_File(surface, bottom_bc, output_file);
        }

        public void Generate_3D_FlexPDE_File(ThreeD_SchrodingerPoissonSolver.Experiment exp, double surface)
        {
            string output_file = "split_gate_3d.pde";

            Check_Output_Filename(ref output_file);

            throw new NotImplementedException();

            //ThreeD_SchrodingerPoissonSolver.ThreeD_PoissonSolver temp_poissolv = new ThreeD_SchrodingerPoissonSolver.ThreeD_PoissonSolver(exp, true, output_file, null, 0.0);
            //temp_poissolv.Create_FlexPDE_File(600, 400, surface, bottom_bc, output_file);
        }

        /// <summary>
        /// checks to see whether the file already exists and asks the user if it should be replaced if it does
        /// </summary>
        void Check_Output_Filename(ref string output_file)
        {
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
        }
    }
}
