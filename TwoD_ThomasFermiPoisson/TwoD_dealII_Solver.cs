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
    class TwoD_dealII_Solver : Potential_Base
    {
        Experiment exp;
        string densdopent_filename;
        string densderiv_filename;
        string pot_filename;
        string new_pot_filename;

        string gphi_filename = "gphi.dat";

        public TwoD_dealII_Solver(Experiment exp, bool using_external_code, string external_input, string external_location, double tol)
            : base(using_external_code, external_input, external_location, tol)
        {
            this.exp = exp;
            this.dens_filename = "car_dens.dat";
            this.densdopent_filename = "dens_2D_dopents.dat";
            this.densderiv_filename = "rho_prime.dat";
            this.pot_filename = "phi.dat";
            this.new_pot_filename = "new_phi.dat";
        }

        protected override string[] Trim_Potential_File(string[] lines)
        {
            // deal.II code should return only potential on grid
            return lines;
        }

        protected override Band_Data Parse_Potential(string[] data)
        {
            string[] new_data = Trim_Potential_File(data);
            return Band_Data.Parse_Band_Data(new_data, exp.Ny_Dens, exp.Nz_Dens);
        }

        protected override Band_Data Get_ChemPot_On_Regular_Grid(Band_Data density)
        {
            throw new NotImplementedException();
        }

        protected override void Save_Density_Data(Band_Data density, string input_file_name)
        {
            density.Save_Data(input_file_name);
        }

        public override Band_Data Calculate_Laplacian(Band_Data input_data)
        {
            DoubleMatrix result = new DoubleMatrix(exp.Ny_Dens, exp.Nz_Dens);
            DoubleMatrix data = input_data.mat;

            for (int i = 1; i < exp.Ny_Dens - 1; i++)
                for (int j = 1; j < exp.Nz_Dens - 1; j++)
                {
                    double pos_y = i * exp.Dy_Dens + exp.Ymin_Dens;
                    double pos_z = j * exp.Dz_Dens + exp.Zmin_Dens;

                    // the factors multiplying the Laplacian in the transverse direction
                    double factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_y + exp.Dy_Dens, pos_z).Permitivity / (exp.Dy_Dens * exp.Dy_Dens);
                    double factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_y - exp.Dy_Dens, pos_z).Permitivity / (exp.Dy_Dens * exp.Dy_Dens);
                    result[i, j] = (factor_minus * data[i - 1, j] + factor_plus * data[i + 1, j] - (factor_plus + factor_minus) * data[i, j]);

                    // and in the growth direction
                    factor_plus = Geom_Tool.GetLayer(exp.Layers, pos_y, pos_z + exp.Dz_Dens).Permitivity / (exp.Dz_Dens * exp.Dz_Dens);
                    factor_minus = Geom_Tool.GetLayer(exp.Layers, pos_y, pos_z - exp.Dz_Dens).Permitivity / (exp.Dz_Dens * exp.Dz_Dens);
                    result[i, j] += (factor_minus * data[i, j - 1] + factor_plus * data[i, j + 1] - (factor_plus + factor_minus) * data[i, j]);
                }

            return new Band_Data(result);
        }

        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data rhs)
        {
            Save_Density_Data(rho_prime.Spin_Summed_Data, densderiv_filename);

            Run_External_Code("x.dat", true);

            string[] lines = File.ReadAllLines("x.dat");
            string[] data = Trim_Potential_File(lines);

            // return chemical potential using mu = - E_c = q_e * phi where E_c is the conduction band edge
            return Physics_Base.q_e * Parse_Potential(data);
        }

        public void Run_External_Code(string filename, bool newton_step)
        {
            Console.WriteLine("Changing state of code to solve Newton step");
            // do something

            Run_External_Code(filename);
        }
    }
}
