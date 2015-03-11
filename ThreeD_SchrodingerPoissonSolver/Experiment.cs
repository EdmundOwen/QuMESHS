using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using CenterSpace.NMath.Core;
using TwoD_ThomasFermiPoisson;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;
using System.Diagnostics;
using System.IO;

namespace ThreeD_SchrodingerPoissonSolver
{
    public class Experiment : Experiment_Base
    {
        double top_V, split_V, surface_charge;
        double split_width, split_length, top_length;
        double t_damp = 0.8, t_min = 1e-3;

        SpinResolved_Data dens_1d;
        ThreeD_ThomasFermiSolver dens_solv;
        ThreeD_PoissonSolver pois_solv;

        public void Initialise_Experiment(Dictionary<string, object> input_dict)
        {
            Console.WriteLine("Initialising Experiment");

            // simulation domain inputs
            Get_From_Dictionary<double>(input_dict, "dx", ref dx_dens); dx_pot = dx_dens;
            Get_From_Dictionary<double>(input_dict, "dy", ref dy_dens); dy_pot = dy_dens;
            Get_From_Dictionary<double>(input_dict, "dz", ref dz_dens); dz_pot = dz_dens;

            Get_From_Dictionary(input_dict, "nx", ref nx_dens); nx_pot = nx_dens - 1;// subtract one from the number of potential points as this is only used to calculate the size of the domain using
            Get_From_Dictionary(input_dict, "ny", ref ny_dens); ny_pot = ny_dens - 1;// lx = nx * dx so on a grid, this is one too many
            Get_From_Dictionary(input_dict, "nz", ref nz_dens); nz_pot = nz_dens;

            // but try to get the specific values
            Get_From_Dictionary<double>(input_dict, "dx_dens", ref dx_dens, true);
            Get_From_Dictionary<double>(input_dict, "dy_dens", ref dy_dens, true);
            Get_From_Dictionary<double>(input_dict, "dz_dens", ref dz_dens, true);
            Get_From_Dictionary<double>(input_dict, "dx_pot", ref dx_pot, true);
            Get_From_Dictionary<double>(input_dict, "dy_pot", ref dy_pot, true);
            Get_From_Dictionary<double>(input_dict, "dz_pot", ref dz_pot, true);

            Get_From_Dictionary(input_dict, "nx_dens", ref nx_dens, true);
            Get_From_Dictionary(input_dict, "ny_dens", ref ny_dens, true);
            Get_From_Dictionary(input_dict, "nz_dens", ref nz_dens, true);
            Get_From_Dictionary(input_dict, "nx_pot", ref nx_pot, true);
            Get_From_Dictionary(input_dict, "ny_pot", ref ny_pot, true);
            Get_From_Dictionary(input_dict, "nz_pot", ref nz_pot, true);

            // physics parameters are done by the base method
            base.Initialise(input_dict);

            // gate voltages
            Get_From_Dictionary<double>(input_dict, "top_V", ref top_V);
            Get_From_Dictionary<double>(input_dict, "split_V", ref split_V);

            // and split gate dimensions
            Get_From_Dictionary<double>(input_dict, "split_width", ref split_width);
            Get_From_Dictionary<double>(input_dict, "split_length", ref split_length);
            Get_From_Dictionary<double>(input_dict, "top_length", ref top_length);

            if (input_dict.ContainsKey("hot_start")) hot_start = (bool)input_dict["hot_start"];
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

                this.carrier_density = new SpinResolved_Data(Band_Data.Parse_Band_Data(spin_up_data, Nx_Dens, Ny_Dens, Nz_Dens), Band_Data.Parse_Band_Data(spin_down_data, Nx_Dens, Ny_Dens, Nz_Dens));

                // and surface charge density
                try { surface_charge = double.Parse(File.ReadAllLines((string)input_dict["surface_charge_file"])[0]); }
                catch (KeyNotFoundException key_e) { throw new Exception("Error - Are the file names for the hot start data included in the input file?\n" + key_e.Message); }
            }
            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            else if (input_dict.ContainsKey("SpinResolved_Density"))
            {
                this.carrier_density = new SpinResolved_Data(new Band_Data(nx_dens, ny_dens, nz_dens, 0.0), new Band_Data(nx_dens, ny_dens, nz_dens, 0.0));
                SpinResolved_Data tmp_charge_1d_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
                dens_1d = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_dens)), new Band_Data(new DoubleVector(nz_dens)));
                int z_offset = (int)Math.Abs((Zmin_Pot - Zmin_Dens) / Dz_Pot);

                for (int i = 0; i < nz_dens; i++)
                {
                    dens_1d.Spin_Up.vec[i] = tmp_charge_1d_density.Spin_Up.vec[z_offset + i];
                    dens_1d.Spin_Down.vec[i] = tmp_charge_1d_density.Spin_Down.vec[z_offset + i];
                }

                // this is the charge density modulation in the (x, y) plane so, initially, just put it in uniformly
                for (int k = 0; k < nz_dens; k++)
                    for (int i = 0; i < nx_dens; i++)
                        for (int j = 0; j < ny_dens; j++)
                        {
                            carrier_density.Spin_Up.vol[k][i, j] = dens_1d.Spin_Up.vec[k];
                            carrier_density.Spin_Down.vol[k][i, j] = dens_1d.Spin_Down.vec[k];
                        }

                Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            }
            else
            {
                this.carrier_density = new SpinResolved_Data(new Band_Data(nx_dens, ny_dens, nz_dens, 0.0), new Band_Data(nx_dens, ny_dens, nz_dens, 0.0));
                Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            }

            // and try to get the dopent density data
            if (input_dict.ContainsKey("Dopent_Density"))
            {
                this.dopent_density = new SpinResolved_Data(new Band_Data(nx_dens, ny_dens, nz_dens, 0.0), new Band_Data(nx_dens, ny_dens, nz_dens, 0.0));
                SpinResolved_Data tmp_1d_dopdens = (SpinResolved_Data)input_dict["Dopent_Density"];

                int offset_min = (int)Math.Abs((Zmin_Pot - Zmin_Dens) / Dz_Pot);
                // this is the charge density modulation in the (x, y) plane so, initially, just put it in uniformly
                for (int k = 0; k < nz_dens; k++)
                    for (int i = 0; i < nx_dens; i++)
                        for (int j = 0; j < ny_dens; j++)
                        {
                            dopent_density.Spin_Up.vol[k][i, j] = tmp_1d_dopdens.Spin_Up.vec[k + offset_min];
                            dopent_density.Spin_Down.vol[k][i, j] = tmp_1d_dopdens.Spin_Down.vec[k + offset_min];
                        }
            }
            else
                this.dopent_density = new SpinResolved_Data(new Band_Data(nx_dens, ny_dens, nz_dens, 0.0), new Band_Data(nx_dens, ny_dens, nz_dens, 0.0));

            // and instantiate their derivatives
            carrier_density_deriv = new SpinResolved_Data(new Band_Data(nx_dens, ny_dens, nz_dens, 0.0), new Band_Data(nx_dens, ny_dens, nz_dens, 0.0));
            dopent_density_deriv = new SpinResolved_Data(new Band_Data(nx_dens, ny_dens, nz_dens, 0.0), new Band_Data(nx_dens, ny_dens, nz_dens, 0.0));

            // calculate the z-position of the maximum of dens_1d and use this as z_2DEG
            double z_2DEG = zmin_dens + dz_dens * (double)Array.IndexOf(dens_1d.Spin_Summed_Data.vec.ToArray(), dens_1d.Spin_Summed_Data.vec.Min());

            // create charge density solver and calculate boundary conditions
            dens_solv = new ThreeD_ThomasFermiSolver(this);

            // initialise potential solver
            bool with_smoothing = false;
            Get_From_Dictionary<bool>(input_dict, "with_smoothing", ref with_smoothing, true);
            pois_solv = new ThreeD_PoissonSolver(this, using_flexPDE, external_input, external_location, tol);
            pois_solv.ZDens = dens_1d;
            pois_solv.Z_2DEG = z_2DEG;

            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            // calculate the bare potential
            Console.WriteLine("Calculating bare potential");
            pois_solv.Set_Boundary_Conditions(top_V, split_V, top_length, split_width, split_length, bottom_V, surface_charge);
            chem_pot = pois_solv.Get_Chemical_Potential(0.0 * carrier_density.Spin_Summed_Data);
            Console.WriteLine("Saving bare potential");
            (Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, dz_dens, nx_dens, ny_dens, nz_dens, xmin_dens, ymin_dens, zmin_dens) - chem_pot).Save_Data("bare_pot.dat");
            Console.WriteLine("Bare potential saved");

            // and then run the DFT solver at the base temperature over a limited range
            //TwoD_DFTSolver dft_solv = new TwoD_DFTSolver(this);
            ThreeD_EffectiveBandSolver dft_solv = new ThreeD_EffectiveBandSolver(this);

            bool converged = false;
            int no_runs = 2500;
            if (no_dft)
                dft_solv.DFT_Mixing_Parameter = 0.0;
            else
                dft_solv.DFT_Mixing_Parameter = 0.1;
            // start without dft if carrier density is empty
            if (carrier_density.Spin_Summed_Data.Min() == 0.0)
                dft_solv.DFT_Mixing_Parameter = 0.0;

            converged = Run_Iteration_Routine(dft_solv, tol, no_runs);

            // initialise output solvers
            ThreeD_ThomasFermiSolver final_dens_solv = new ThreeD_ThomasFermiSolver(this);
            ThreeD_PoissonSolver final_pois_solv = new ThreeD_PoissonSolver(this, using_flexPDE, external_input, external_location, tol);

            // save final density out
            carrier_density.Spin_Summed_Data.Save_Data("dens_3D.dat");
            carrier_density.Spin_Up.Save_Data("dens_3D_up.dat");
            carrier_density.Spin_Down.Save_Data("dens_3D_down.dat");

            // save surface charge
            StreamWriter sw = new StreamWriter("surface_charge.dat"); sw.WriteLine(surface_charge.ToString()); sw.Close();
            // save eigen-energies
            DoubleVector energies = dft_solv.Get_EnergyLevels(layers, chem_pot);
            StreamWriter sw_e = new StreamWriter("energies.dat");
            for (int i = 0; i < energies.Length; i++)
                sw_e.WriteLine(energies[i]);
            sw_e.Close();

            final_dens_solv.Output(carrier_density, "carrier_density.dat");
            final_dens_solv.Output(carrier_density - dft_solv.Get_ChargeDensity(layers, carrier_density, dopent_density, chem_pot), "density_error.dat");
            (Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, dz_dens, nx_dens, ny_dens, nz_dens, xmin_dens, ymin_dens, zmin_dens) - chem_pot).Save_Data("potential.dat");
            Band_Data pot_exc = dft_solv.DFT_diff(carrier_density) + Physics_Base.Get_XC_Potential(carrier_density);
            pot_exc.Save_Data("xc_pot.dat");
            (Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, dz_dens, nx_dens, ny_dens, nz_dens, xmin_dens, ymin_dens, zmin_dens) - chem_pot + pot_exc).Save_Data("pot_KS.dat");
            Band_Data ks_ke = dft_solv.Get_KS_KE(layers, chem_pot);
            ks_ke.Save_Data("ks_ke.dat");

            throw new NotImplementedException();
        }

        bool Run_Iteration_Routine(IDensity_Solve dens_solv, double pot_lim)
        {
            return Run_Iteration_Routine(dens_solv, pot_lim, int.MaxValue);
        }

        double dens_diff_lim = 0.1; // the maximum percentage change in the density required for update of V_xc
        double max_vxc_diff = double.MaxValue; // maximum difference for dft potential... if this increases, the dft mixing parameter is reduced
        double min_dens_diff = 0.02; // minimum bound for the required, percentage density difference for updating the dft potential
        double min_vxc_diff = 0.1; // minimum difference in the dft potential for convergence
        double min_alpha = 0.03; // minimum possible value of the dft mixing parameter
        bool Run_Iteration_Routine(IDensity_Solve dens_solv, double pot_lim, int max_count)
        {
            // calculate initial potential with the given charge distribution
            Console.WriteLine("Calculating initial potential grid");
            pois_solv.Set_Boundary_Conditions(top_V, split_V, top_length, split_width, split_length, bottom_V, surface_charge);
            chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
            Console.WriteLine("Initial grid complete");
            dens_solv.Set_DFT_Potential(carrier_density);
            dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
            dens_solv.Set_DFT_Potential(carrier_density); 

            int count = 0;
            bool converged = false;
            if (!no_dft)
                dens_solv.DFT_Mixing_Parameter = 0.1;
            dens_diff_lim = 0.12;
            while (!converged)
            {
                Stopwatch stpwch = new Stopwatch();
                stpwch.Start();

                // save old density data
                Band_Data dens_old = carrier_density.Spin_Summed_Data.DeepenThisCopy();

                // Get charge rho(phi) (not dopents as these are included as a flexPDE input)
                dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                Set_Edges(carrier_density);

                // Generate an approximate charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi) using the Thomas-Fermi semi-classical method
                SpinResolved_Data rho_prime = dens_solv.Get_ChargeDensity_Deriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);
                Set_Edges(rho_prime);

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                Band_Data gphi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data;
                Set_Edges(gphi);
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, gphi, carrier_density, dens_solv.DFT_diff(carrier_density));
                chem_pot = pois_solv.Chemical_Potential;

                // Calculate optimal damping parameter, t, (but damped damping....)
                if (t == 0.0)
                    t = t_min;

                t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot, x, carrier_density, dopent_density, pois_solv, dens_solv, t_min);
                if (t < 0.0)
                {
                    Console.WriteLine("Iterator has stalled, setting t = 0");
                    t = 0.0;
                }

                // Check convergence
                double[] diff = new double[nx_dens * ny_dens * nz_dens];
                for (int i = 0; i < nx_dens * ny_dens * nz_dens; i++)
                    diff[i] = Math.Abs(gphi[i]);
                double convergence = diff.Sum();

                // and check convergence of density
                Band_Data dens_diff = carrier_density.Spin_Summed_Data - dens_old;
                Band_Data car_dens_spin_summed = carrier_density.Spin_Summed_Data;
                double carrier_dens_min = Math.Abs(car_dens_spin_summed.Min());
                // using the relative absolute density difference
                for (int i = 0; i < dens_diff.Length; i++)
                    // only calculate density difference for densities more than 1% of the maximum value
                    if (Math.Abs(car_dens_spin_summed[i]) > 0.01 * carrier_dens_min)
                        dens_diff[i] = Math.Abs(dens_diff[i] / car_dens_spin_summed[i]);
                    else
                        dens_diff[i] = 0.0;

                // only renew DFT potential when the difference in density has converged and the iterator has done at least 3 iterations
                if (dens_diff.Max() < dens_diff_lim && t > 10.0 * t_min && count > 3)
                {
                    // and set the DFT potential
                    if (dens_solv.DFT_Mixing_Parameter != 0.0)
                        dens_solv.Print_DFT_diff(carrier_density);
                    dens_solv.Set_DFT_Potential(carrier_density);

                    // also... if the difference in the old and new dft potentials is greater than for the previous V_xc update, reduce the dft mixing parameter
                    double current_vxc_diff = Math.Max(dens_solv.DFT_diff(carrier_density).Max(), (-1.0 * dens_solv.DFT_diff(carrier_density).Min()));
                    if (current_vxc_diff > max_vxc_diff && dens_diff_lim / 2.0 > min_dens_diff)
                    {
                        dens_diff_lim /= 2.0;
                        Console.WriteLine("Minimum percentage density difference reduced to " + dens_diff_lim.ToString());
                    }
                    max_vxc_diff = current_vxc_diff;

                    // solution is converged if the density accuracy is better than half the minimum possible value for changing the dft potential
                    if (dens_diff.Max() < min_dens_diff / 2.0 && current_vxc_diff < min_vxc_diff)
                        converged = true;
                }

                // update band energy phi_new = phi_old + t * x
                chem_pot = chem_pot + t * x;
                pois_solv.T = t;

                stpwch.Stop();
                Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString("F") + "\tt = " + t.ToString() + "\ttime = " + stpwch.Elapsed.TotalMinutes.ToString("F"));
                count++;

                // reset the potential if the added potential t * x is too small
                if (converged || count > max_count)
                {
                    File.Copy("split_gate.pg6", "split_gate_final.pg6", true);
                    Console.WriteLine("Maximum potential change at end of iteration was " + Math.Max(t * x.Max(), (-t * x).Max()).ToString());
                    break;
                }
            }

            Console.WriteLine("Iteration complete");
            return converged;
        }

        /// <summary>
        /// sets the edges to be the same as their nearest neighbour.
        /// traditionally, would expect the edges to have zero charge but this will be different for simulations in the plane of the 2deg
        /// </summary>
        /// <param name="carrier_density"></param>
        private void Set_Edges(SpinResolved_Data carrier_density)
        {
            Set_Edges(carrier_density.Spin_Up);
            Set_Edges(carrier_density.Spin_Down);
        }

        private void Set_Edges(Band_Data data)
        {
            for (int k = 0; k < nz_dens; k++)
                for (int i = 0; i < nx_dens; i++)
                {
                    data.vol[k][i, 0] = data.vol[k][i, 1];
                    data.vol[k][i, ny_dens - 1] = data.vol[k][i, ny_dens - 2];
                }

            for (int k = 0; k < nz_dens; k++)
                for (int j = 0; j < ny_dens; j++)
                {
                    data.vol[k][0, j] = data.vol[k][1, j];
                    data.vol[k][nx_dens - 1, j] = data.vol[k][nx_dens - 2, j];
                }

            for (int i = 0; i < nx_dens; i++)
                for (int j = 0; j < ny_dens; j++)
                {

                    data.vol[0][i, j] = data.vol[1][i, j];
                    data.vol[nz_dens - 1][i, j] = data.vol[nz_dens - 2][i, j];
                }
        }
    }
}
