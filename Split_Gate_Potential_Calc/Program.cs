using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using OneD_ThomasFermiPoisson;
using Solver_Bases.Geometry;

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
            input.Add("TF_only", true);             // only to Thomas-Fermi calculations for finding split gate potentials
            input["use_FlexPDE"] = "false";
            Experiment exp = new Experiment();
            exp.Initialise(input);
            exp.Run();

            // convert 1d density to 2d
            double dy = (double)input["dy"]; double dz = (double)input["dz"];
            int ny = (int)(double)input["ny"]; int nz = (int)(double)input["nz"];
            double ymin =  -0.5 * ny * dy;
            double zmin = Solver_Bases.Geometry.Geom_Tool.Get_Zmin(layers);
            input.Add("ymin", ymin); input.Add("zmin", zmin);
            SpinResolved_Data charge_dens_1d = exp.Charge_Density;
            Band_Data charge_dens_2d = Input_Band_Structure.Expand_BandStructure(charge_dens_1d.Spin_Summed_Data.vec, ny);
            charge_dens_2d.Save_2D_Data("dens_2D.dat", dy, dz, ymin, zmin);

            // create a poisson solver and get the possible surface charge
            OneD_PoissonSolver tmp_pois_solv = new OneD_PoissonSolver(dz, nz, layers, false, "", "", 0.0);
            double surface_charge = tmp_pois_solv.Get_Surface_Charge(charge_dens_1d.Spin_Summed_Data, layers);
            Band_Data surface_dens_2d = new Band_Data(new CenterSpace.NMath.Core.DoubleMatrix(ny, nz));
            int zsurface_index = (int)Math.Floor((layers[Geom_Tool.Find_Layer_Above_Surface(layers)].Zmin - Geom_Tool.Get_Zmin(layers)) / dz);
            for (int i = 0; i < ny; i++)
                surface_dens_2d.mat[i, zsurface_index] = surface_charge;
            surface_dens_2d.Save_2D_Data("dens_surface.dat", dy, dz, ymin, zmin);

            // delete unnecessary files
            if (File.Exists("dens_1D.dat"))
                File.Delete("dens_1D.dat");
            if (File.Exists("charge_density.dat"))
                File.Delete("charge_density.dat");
            if (File.Exists("potential.dat"))
                File.Delete("potential.dat");

            // create split gate FlexPDE file (for the moment, only in one dimension)
            SG_File_Generator sg_file_gen = new SG_File_Generator(input);
            sg_file_gen.Generate_2D_FlexPDE_File();
        }
    }
}
