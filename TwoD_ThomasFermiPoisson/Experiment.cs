using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    public class Experiment : Experiment_Base
    {
        double top_V, split_V, surface_charge;

        TwoD_ThomasFermiSolver dens_solv;
        TwoD_PoissonSolver pois_solv;

        public void Initialise_Experiment(Dictionary<string, object> input_dict)
        {
            Console.WriteLine("Initialising Experiment");
            
            // simulation domain inputs
            Get_From_Dictionary<double>(input_dict, "dy", ref dy_dens); dy_pot = dy_dens;
            Get_From_Dictionary<double>(input_dict, "dz", ref dz_dens); dz_pot = dz_dens;

            Get_From_Dictionary(input_dict, "ny", ref ny_dens); ny_pot = ny_dens;
            Get_From_Dictionary(input_dict, "nz", ref nz_dens); nz_pot = nz_dens;

            // but try to get the specific values
            Get_From_Dictionary<double>(input_dict, "dy_dens", ref dy_dens, true);
            Get_From_Dictionary<double>(input_dict, "dz_dens", ref dz_dens, true);
            Get_From_Dictionary<double>(input_dict, "dy_pot", ref dy_pot, true);
            Get_From_Dictionary<double>(input_dict, "dz_pot", ref dz_pot, true);

            Get_From_Dictionary(input_dict, "ny_dens", ref ny_dens, true);
            Get_From_Dictionary(input_dict, "nz_dens", ref nz_dens, true);
            Get_From_Dictionary(input_dict, "ny_pot", ref ny_pot, true);
            Get_From_Dictionary(input_dict, "nz_pot", ref nz_pot, true);

            // physics parameters are done by the base method
            base.Initialise(input_dict);
            
            // gate voltages
            Get_From_Dictionary<double>(input_dict, "top_V", ref top_V);
            Get_From_Dictionary<double>(input_dict, "split_V", ref split_V);
            Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density"))
            {
                this.charge_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));
                SpinResolved_Data tmp_1d_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];

                int offset_min = (int)Math.Round((zmin_dens - zmin_pot) / dz_pot);

                for (int i = 0; i  < ny_dens; i++)
                    for (int j = 0; j < nz_dens; j++)
                    {
                        charge_density.Spin_Up.mat[i, j] = tmp_1d_density.Spin_Up.vec[offset_min + j];
                        charge_density.Spin_Down.mat[i, j] = tmp_1d_density.Spin_Down.vec[offset_min + j];
                    }
            }
            else
                this.charge_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));

            if (input_dict.ContainsKey("dft")) this.TF_only = !bool.Parse((string)input_dict["dft"]);
            if (input_dict.ContainsKey("TF_only")) this.TF_only = bool.Parse((string)input_dict["TF_only"]);

            // create charge density solver and calculate boundary conditions
            dens_solv = new TwoD_ThomasFermiSolver(temperature, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens);
            double bottom_bc = dens_solv.Get_Chemical_Potential(0.0, zmin_pot, layers);

            // initialise potential solver
            pois_solv = new TwoD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            pois_solv.Set_Boundary_Conditions(top_V, split_V, bottom_bc, surface_charge);

            // initialise the band energy as the solution from the input density
            //if (input_dict.ContainsKey("Band_Offset"))
            //{ Band_Data tmp_1d_bandoffset = (Band_Data)input_dict["Band_Offset"]; this.band_offset = Input_Band_Structure.Expand_BandStructure(tmp_1d_bandoffset.vec, ny_dens); }
            //else 
                //band_offset = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data); 

            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            int count = 0;
            while (!dens_solv.Converged)
            {
                Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + temperature.ToString() + "\tConvergence factor: " + dens_solv.Convergence_Factor.ToString());

                // solve the band energy for the given charge  density and mix in with the old band energy
                band_offset = -1.0 * pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data);

                // find the density for this new band offset and blend
                SpinResolved_Data new_density = dens_solv.Get_ChargeDensity(layers, charge_density, band_offset);
                dens_solv.Blend(ref charge_density, new_density, alpha, tol);

                //pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                count++;
            }

            // initialise output solvers
            TwoD_ThomasFermiSolver final_dens_solv = new TwoD_ThomasFermiSolver(temperature, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens);
            TwoD_PoissonSolver final_pois_solv = new TwoD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            charge_density.Spin_Summed_Data.Save_2D_Data("dens_2D.dat", dy_dens, dz_dens, ymin_dens, zmin_dens);

            final_dens_solv.Output(charge_density, "charge_density.dat");
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens) - band_offset, "potential.dat");
        }
    }
}
