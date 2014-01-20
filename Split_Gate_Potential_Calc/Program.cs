using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;
using OneD_ThomasFermiPoisson;

namespace Split_Gate_Potential_Calc
{
    class Program
    {
        static void Main(string[] args)
        {
            // Input simulation data
            string input_file;
            if (args.Length == 0)
            {
                Console.WriteLine("Input file name:");
                input_file = Console.ReadLine();
            }
            else
                input_file = args[0];

            if (!File.Exists(input_file))
                throw new FileNotFoundException("Error - input file \"" + input_file + "\" not found!");

            Dictionary<string, object> input = new Dictionary<string,object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref input, input_file);

            // Input band structure data
            string band_file;
            if (args.Length != 2)
            {
                Console.WriteLine("Band structure file name:");
                band_file = Console.ReadLine();
            }
            else
                band_file = args[1];

            if (!File.Exists(band_file))
                throw new FileNotFoundException("Error - input file \"" + band_file + "\" containing band structure information not found!");

            ILayer[] layers = Input_Band_Structure.Get_Layers(band_file);
            input.Add("Layers", layers);

            // calculate 1D potentials
            input.Add("TF_only", "true");             // only to Thomas-Fermi calculations for finding split gate potentials
            input["use_FlexPDE"] = "false";
            Experiment exp = new Experiment();
            exp.Initialise(input);
            exp.Run();

            // create a poisson solver and get the possible surface charge
            double dz = (double)input["dz"]; int nz = (int)(double)input["nz"];
            double zmin = Solver_Bases.Geometry.Geom_Tool.Get_Zmin(layers);
            Band_Data band_offset = exp.Band_Offset;
            OneD_PoissonSolver tmp_pois_solv = new OneD_PoissonSolver(exp, false, "", "", 0.0);
            double surface_charge = tmp_pois_solv.Get_Surface_Charge(band_offset, layers);

            // save the 1D data
            SpinResolved_Data charge_dens_1d = exp.Charge_Density;
            charge_dens_1d.Spin_Summed_Data.Save_1D_Data("dens_1D.dat", dz, zmin);
            input.Add("oned_dens_data", charge_dens_1d.Spin_Summed_Data.vec);
            
            // delete unnecessary files
            //if (File.Exists("dens_1D.dat"))
            //    File.Delete("dens_1D.dat");
            if (File.Exists("charge_density.dat"))
                File.Delete("charge_density.dat");
            if (File.Exists("potential.dat"))
                File.Delete("potential.dat");

            // add parameters to input dictionary
            double dy = (double)input["dy"]; 
            int ny = (int)(double)input["ny"]; 
            double ymin = -0.5 * ny * dy;
            input.Add("ymin", ymin); input.Add("zmin", zmin);

            SG_File_Generator sg_file_gen = new SG_File_Generator(input);

            // create relevant flexpde file
            if (bool.Parse((string)input["is_2D"]))
            {
                Band_Data charge_dens_2d = Input_Band_Structure.Expand_BandStructure(charge_dens_1d.Spin_Summed_Data.vec, ny);
                charge_dens_2d.Save_2D_Data("dens_2D.dat", dy, dz, ymin, zmin);

                // create split gate FlexPDE file (for the moment, only in one dimension)
                TwoD_ThomasFermiPoisson.Experiment exp_2d = new TwoD_ThomasFermiPoisson.Experiment();
                exp_2d.Initialise_Experiment(input);

                sg_file_gen.Generate_2D_FlexPDE_File(exp_2d, surface_charge * dz);
            }
            else
            {
                double dx = (double)input["dx"]; int nx = (int)(double)input["nx"];
                double xmin = -0.5 * nx * dx;
                input.Add("xmin", xmin);

                Band_Data charge_dens_3d = Input_Band_Structure.Expand_BandStructure(charge_dens_1d.Spin_Summed_Data.vec, nx, ny);
                charge_dens_3d.Save_3D_Data("dens_3D.dat", dx, dy, dz, xmin, ymin, zmin);

                // create split gate FlexPDE file (for the moment, only in one dimension)
                sg_file_gen.Generate_3D_FlexPDE_File(surface_charge * dz);
            }
        }
    }
}
