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
        string dens_filename = "car_dens.dat";
        string densderiv_filename = "rho_prime.dat";
        string gphi_filename = "gphi.dat";
                
        Band_Data chempot;

        public TwoD_dealII_Solver(Experiment exp, bool using_external_code, Dictionary<string, object> input)
            : base(using_external_code)
        {
            this.exp = exp;

            if (!input.ContainsKey("initcalc_location") || !input.ContainsKey("newton_location"))
                throw new Exception("Error - To use deal.II you must provide the location for calculating the initial potential in \"initcalc_location\" and the newton step calculation in \"newton_location\"");

            this.initcalc_location = (string)input["initcalc_location"];
            this.newton_location = (string)input["newton_location"];
        }

        public override void Initiate_Poisson_Solver(Dictionary<string, double> device_dimensions, Dictionary<string, double> boundary_conditions)
        {

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

        protected override void Save_Data(Band_Data density, string input_file_name)
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

        protected override Band_Data Get_ChemPot_From_External(Band_Data density)
        {
            Save_Data(density, dens_filename);

            chempot = Get_Data_From_External(initcalc_location, initcalc_result_filename);
            return chempot;
        }

        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi)
        {
            Save_Data(rho_prime.Spin_Summed_Data, densderiv_filename);
            Save_Data(gphi, gphi_filename);

            Band_Data x = Get_Data_From_External(newton_location, newton_result_filename);
            chempot += base.T * x;
            return x;
        }

        public override Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data car_dens, Band_Data dft_diff)
        {
            Save_Data(dft_diff + Physics_Base.Get_XC_Potential(car_dens), "xc_pot_calc.dat");
            Save_Data(car_dens.Spin_Summed_Data, dens_filename);
            Save_Data(Physics_Base.Get_XC_Potential(car_dens), xc_pot_filename);

            return Calculate_Newton_Step(rho_prime, gphi);
        }

        public override Band_Data Chemical_Potential
        {
            get { return chempot; }
        }
    }
}
