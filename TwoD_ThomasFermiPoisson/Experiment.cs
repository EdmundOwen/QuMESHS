using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;
using Solver_Bases.Mixing_Schedulers;

namespace TwoD_ThomasFermiPoisson
{
    public class Experiment : Experiment_Base
    {
        double top_V, split_V, surface_charge;
        double split_width;

        TwoD_ThomasFermiSolver dens_solv;
        TwoD_PoissonSolver pois_solv;
        IScheduler scheduler;

        double max_dens_init;
        Band_Data[] pois_results;
        Band_Data pois_background;

        public void Initialise_Experiment(Dictionary<string, object> input_dict)
        {
            Console.WriteLine("Initialising Experiment");

            // mixing parameter scheduler
            if (!input_dict.ContainsKey("scheduler_file"))
                input_dict.Add("scheduler_file", "none");
            scheduler = Mixing_Scheduler_Factory.Get_Scheduler((string)input_dict["scheduler_file"]);
            
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

            // and split gate dimensions
            Get_From_Dictionary<double>(input_dict, "split_width", ref split_width);

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            bool hot_start = false;
            if (input_dict.ContainsKey("hot_start")) hot_start = bool.Parse((string)input_dict["hot_start"]);
            if (hot_start)
            {
                // load (spin-resolved) density data
                string[] spin_up_data, spin_down_data;
                try
                {
                    spin_up_data = File.ReadAllLines((string)input_dict["spin_up_file"]);
                    spin_down_data = File.ReadAllLines((string)input_dict["spin_down_file"]);
                }
                catch (KeyNotFoundException key_e)
                { throw new Exception("Error - Are the file names for the hot start data included in the input file?\n" + key_e.Message); }

                this.carrier_density = new SpinResolved_Data(Band_Data.Parse_Band_Data(spin_up_data, Ny_Dens, Nz_Dens), Band_Data.Parse_Band_Data(spin_down_data, Ny_Dens, Nz_Dens));

                // and surface charge density
                try { surface_charge = double.Parse(File.ReadAllLines((string)input_dict["surface_charge_file"])[0]); }
                catch (KeyNotFoundException key_e) {throw new Exception("Error - Are the file names for the hot start data included in the input file?\n" + key_e.Message); }
            }
            else if (input_dict.ContainsKey("SpinResolved_Density"))
            {
                this.carrier_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));

                SpinResolved_Data tmp_1d_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];

                int offset_min = (int)Math.Round((zmin_dens - zmin_pot) / dz_pot);

                for (int i = 0; i < ny_dens; i++)
                    for (int j = 0; j < nz_dens; j++)
                    {
                        // do not add anything to the density if on the edge of the domain
                        if (i == 0 || i == ny_dens - 1 || j == 0 || j == nz_dens - 1)
                            continue;

                        // linearly interpolate
                        int j_init = (int)Math.Floor(j * Dz_Dens / Dz_Pot);
                        carrier_density.Spin_Up.mat[i, j] = tmp_1d_density.Spin_Up.vec[offset_min + j_init] + (tmp_1d_density.Spin_Up.vec[offset_min + j_init + 1] - tmp_1d_density.Spin_Up.vec[offset_min + j_init]) * Dz_Dens / Dz_Pot;
                        carrier_density.Spin_Down.mat[i, j] = tmp_1d_density.Spin_Down.vec[offset_min + j_init] + (tmp_1d_density.Spin_Down.vec[offset_min + j_init + 1] - tmp_1d_density.Spin_Down.vec[offset_min + j_init]) * Dz_Dens / Dz_Pot;
                    }

                Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            }
            else
            {
                this.carrier_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));
                Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            }
            
            if (input_dict.ContainsKey("dft")) this.TF_only = !bool.Parse((string)input_dict["dft"]);
            if (input_dict.ContainsKey("TF_only")) this.TF_only = bool.Parse((string)input_dict["TF_only"]);

            // create charge density solver and calculate boundary conditions
            dens_solv = new TwoD_ThomasFermiSolver(temperature, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens);
            double bottom_V = dens_solv.Get_Chemical_Potential(0.0, zmin_pot, layers) / (Physics_Base.q_e * Physics_Base.energy_V_to_meVpzC);

            // initialise potential solver
            pois_solv = new TwoD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);

            // calculate what the maximum magnitude density is
            max_dens_init = Math.Max(Math.Abs(carrier_density.Spin_Summed_Data.mat.Min()), carrier_density.Spin_Summed_Data.mat.Max());

            Console.WriteLine("Calculating potentials");
            pois_results = new Band_Data[ny_dens * nz_dens];
            // calculate what the potential is with no charges
            Band_Data pois_nocharge = pois_solv.Calculate_Point_Potential(top_V, split_V, split_width, bottom_V, surface_charge, 1, 1, 0.0);

            // and calculate the background by setting the point charge magnitude to zero
            pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);
            pois_background = pois_solv.Get_Chemical_Potential(0.0 * carrier_density.Spin_Summed_Data);

            // and use this value as the magnitude of a delta function for the potential for each of the density points minus the external potential
            for (int i = 0; i < ny_dens; i++)
                for (int j = 0; j < nz_dens; j++)
                {
                    Band_Data tmp = pois_solv.Calculate_Point_Potential(top_V, split_V, split_width, bottom_V, surface_charge, i, j, max_dens_init);
                    if (tmp != null)
                        pois_results[i * nz_dens + j] = tmp - pois_nocharge;
                    else
                        pois_results[i * nz_dens + j] = null;
                    Console.WriteLine("Potential " + (i * nz_dens + j).ToString() + "/" + (ny_dens * nz_dens).ToString() + " calculated");
                }

            Console.WriteLine("Potentials calculated");

            // initialise the band energy as the solution from the input density
            //if (input_dict.ContainsKey("Band_Offset"))
            //{ Band_Data tmp_1d_bandoffset = (Band_Data)input_dict["Band_Offset"]; this.band_offset = Input_Band_Structure.Expand_BandStructure(tmp_1d_bandoffset.vec, ny_dens); }
            //else 
                //band_offset = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data); 

            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            SpinResolved_Data old_carrier_density = carrier_density.DeepenThisCopy();

            int count = 0;

            if (TF_only)
            {
                while (!dens_solv.Converged || count < 25)
                {
                    Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + temperature.ToString() + "\tConvergence factor: " + dens_solv.Convergence_Factor.ToString());

                    // solve the chemical potential for the given charge  density and mix in with the old chemical potential
                    chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);

                    // find the density for this new chemical potential and blend
                    SpinResolved_Data new_carrier_density = dens_solv.Get_ChargeDensity(layers, carrier_density, chem_pot);
                    alpha = scheduler.Get_Mixing_Parameter(count, "alpha"); double zeta = scheduler.Get_Mixing_Parameter(count, "zeta");
                    dens_solv.Blend(ref carrier_density, ref old_carrier_density, new_carrier_density, alpha, zeta, tol);

                    count++;
                }
            }

            // and then run the DFT solver at the base temperature over a limited range
            TwoD_DFTSolver dft_solv = new TwoD_DFTSolver(this);
            //TwoD_EffectiveBandSolver dft_solv = new TwoD_EffectiveBandSolver(this);

            count = 0;
            while ((!dft_solv.Converged || count < 25) && TF_only != true)
            {
                Console.WriteLine("Iteration: " + count.ToString() + "\tTemp: " + temperature.ToString() + "\tConvergence factor: " + dft_solv.Convergence_Factor.ToString());

                // solve the chemical potential for the given charge  density
                chem_pot = pois_solv.Get_Chemical_Potential(pois_background, pois_results, carrier_density.Spin_Summed_Data, 1.0 / max_dens_init);

                // find the density for this new chemical potential and blend
                Band_Data dft_chem_pot = Get_Potential(ref chem_pot, layers);
                SpinResolved_Data new_carrier_density = dft_solv.Get_ChargeDensity(layers, carrier_density, dft_chem_pot);
                alpha = scheduler.Get_Mixing_Parameter(count, "alpha"); double zeta = scheduler.Get_Mixing_Parameter(count, "zeta");
                dft_solv.Blend(ref carrier_density, ref old_carrier_density, new_carrier_density, alpha, zeta, tol);
                
                count++;
            }

            // initialise output solvers
            TwoD_ThomasFermiSolver final_dens_solv = new TwoD_ThomasFermiSolver(temperature, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens);
            TwoD_PoissonSolver final_pois_solv = new TwoD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            carrier_density.Spin_Summed_Data.Save_Data("dens_2D_raw.dat");
            carrier_density.Spin_Up.Save_Data("dens_2D_up_raw.dat");
            carrier_density.Spin_Down.Save_Data("dens_2D_down_raw.dat");

            // save surface charge
            StreamWriter sw = new StreamWriter("surface_charge.dat"); sw.WriteLine(surface_charge.ToString()); sw.Close();
            // save eigen-energies
            DoubleVector energies = dft_solv.Get_EnergyLevels(layers, carrier_density, Get_Potential(ref chem_pot, layers));
            StreamWriter sw_e = new StreamWriter("energies.dat");
            for (int i = 0; i < energies.Length; i++)
                sw_e.WriteLine(energies[i]);
            sw_e.Close();

            final_dens_solv.Output(carrier_density, "carrier_density.dat");
            final_dens_solv.Output(carrier_density - dft_solv.Get_ChargeDensity(layers, carrier_density, Get_Potential(ref chem_pot, layers)), "density_error.dat");
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens) - chem_pot, "potential.dat");
        }

        Band_Data Get_Potential(ref Band_Data chem_pot, ILayer[] layers)
        {
            Band_Data result = new Band_Data(new DoubleMatrix(ny_dens, nz_dens));

            for (int i = 0; i < ny_dens; i++)
                for (int j = 0; j < nz_dens; j++)
                {
                    double pos_y = ymin_dens + i * dy_dens;
                    double pos_z = zmin_dens + j * dz_dens;

                    double band_gap = Geom_Tool.GetLayer(layers, pos_y, pos_z).Band_Gap;
                    result.mat[i, j] = 0.5 * band_gap - chem_pot.mat[i, j];
                }

            return result;
        }
    }
}
