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
        double t_damp = 0.8, t_min = 1e-3;

        IPoisson_Solve pois_solv;

        public override void Initialise(Dictionary<string, object> input_dict)
        {
            Console.WriteLine("Initialising Experiment");

            // simulation domain inputs
            Get_From_Dictionary<double>(input_dict, "dx", ref dx_dens); dx_pot = dx_dens;
            Get_From_Dictionary<double>(input_dict, "dy", ref dy_dens); dy_pot = dy_dens;
            Get_From_Dictionary<double>(input_dict, "dz", ref dz_dens); dz_pot = dz_dens;

            Get_From_Dictionary(input_dict, "nx", ref nx_dens); nx_pot = nx_dens;
            Get_From_Dictionary(input_dict, "ny", ref ny_dens); ny_pot = ny_dens;
            Get_From_Dictionary(input_dict, "nz", ref nz_dens); nz_pot = nz_dens;

            // physics parameters are done by the base method (the base will also try to get specific parameters detailed in the input files)
            base.Initialise(input_dict);

            // and initialise the data classes for density, its derivatives and the chemical potential
            Initialise_DataClasses(input_dict);

            // initialise potential solver
            if (using_flexPDE)
                pois_solv = new ThreeD_PoissonSolver(this, using_flexPDE, input_dict);
            else if (using_dealii)
                pois_solv = new ThreeD_dealII_Solver(this, using_dealii, input_dict);
            else
                throw new NotImplementedException("Error - Must use either FlexPDE or deal.II for 2D potential solver!");

            device_dimensions.Add("interface_depth", Layers[1].Zmax);
            pois_solv.Initiate_Poisson_Solver(device_dimensions, boundary_conditions);

            // and load the output suffix for identification of output files
            Get_From_Dictionary<string>(input_dict, "output_suffix", ref output_suffix, true);

            Console.WriteLine("Experimental parameters initialised");
        }

        protected override void Initialise_DataClasses(Dictionary<string, object> input_dict)
        {
            // initialise data classes for the density and chemical potential
            this.carrier_density = new SpinResolved_Data(nx_dens, ny_dens, nz_dens);
            this.dopent_density = new SpinResolved_Data(nx_dens, ny_dens, nz_dens);
            this.chem_pot = new Band_Data(nx_dens, ny_dens, nz_dens, 0.0);

            if (hot_start)
            {
                Initialise_from_Hot_Start(input_dict);
            }
            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            else if (initialise_with_1D_data)
            {
                Initialise_from_1D(input_dict);
            }
            else if (initialise_from_restart)
            {
                Initialise_from_Restart(input_dict);
            }
            else
            {
                boundary_conditions.Add("surface", (double)input_dict["surface_charge"]);
            }

            // and instantiate their derivatives
            carrier_density_deriv = new SpinResolved_Data(nx_dens, ny_dens, nz_dens);
            dopent_density_deriv = new SpinResolved_Data(nx_dens, ny_dens, nz_dens);
        }
        
        public override void Run()
        {
            if (!initialise_from_restart)
            {
                // calculate the bare potential
                Console.WriteLine("Calculating bare potential");
                chem_pot = Physics_Base.q_e * pois_solv.Get_Potential(0.0 * carrier_density.Spin_Summed_Data);
                Console.WriteLine("Saving bare potential");
                (Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, dz_dens, nx_dens, ny_dens, nz_dens, xmin_dens, ymin_dens, zmin_dens) - chem_pot).Save_Data("bare_pot.dat");
                Console.WriteLine("Bare potential saved");

                // if the initial carrier density was not zero, recalculate the chemical potential
                if (carrier_density.Spin_Summed_Data.Max() != 0.0 || carrier_density.Spin_Summed_Data.Min() != 0.0)
                    chem_pot = Physics_Base.q_e * pois_solv.Get_Potential(carrier_density.Spin_Summed_Data);
            }

            ThreeD_ThomasFermiSolver dft_solv = new ThreeD_ThomasFermiSolver(this);
          //  ThreeD_EffectiveBandSolver dft_solv = new ThreeD_EffectiveBandSolver(this);
          //  TwoplusOneD_ThomasFermiSolver dft_solv = new TwoplusOneD_ThomasFermiSolver(this);

            bool converged = false;
            // start without dft if carrier density is empty
            if (no_dft || carrier_density.Spin_Summed_Data.Min() == 0.0)
                dft_solv.DFT_Mixing_Parameter = 0.0;
            else
                dft_solv.DFT_Mixing_Parameter = 0.1;

            converged = Run_Iteration_Routine(dft_solv, pois_solv, tol, max_iterations);

            // save surface charge
            StreamWriter sw = new StreamWriter("surface_charge.dat"); sw.WriteLine(boundary_conditions["surface"].ToString()); sw.Close();
            // save eigen-energies
            DoubleVector energies = dft_solv.Get_EnergyLevels(layers, chem_pot);
            StreamWriter sw_e = new StreamWriter("energies.dat");
            for (int i = 0; i < energies.Length; i++)
                sw_e.WriteLine(energies[i]);
            sw_e.Close();

            dft_solv.Output(carrier_density, "carrier_density.dat");
            dft_solv.Output(carrier_density - dft_solv.Get_ChargeDensity(layers, carrier_density, dopent_density, chem_pot), "density_error.dat");
            (Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, dz_dens, nx_dens, ny_dens, nz_dens, xmin_dens, ymin_dens, zmin_dens) - chem_pot).Save_Data("potential.dat");
            Band_Data pot_exc = dft_solv.DFT_diff(carrier_density) + Physics_Base.Get_XC_Potential(carrier_density);
            pot_exc.Save_Data("xc_pot.dat");
            (Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, dz_dens, nx_dens, ny_dens, nz_dens, xmin_dens, ymin_dens, zmin_dens) - chem_pot + pot_exc).Save_Data("pot_KS.dat");
//            Band_Data ks_ke = dft_solv.Get_KS_KE(layers, chem_pot);
//            ks_ke.Save_Data("ks_ke.dat");

            // clean up intermediate data files
            File.Delete("phi.dat");
            File.Delete("new_phi.dat");
            File.Delete("x.dat");
            File.Delete("y.dat");
            File.Delete("gphi.dat");
            File.Delete("car_dens.dat");
            File.Delete("rho_prime.dat");
            File.Delete("xc_pot.dat");
            File.Delete("xc_pot_calc.dat");
            File.Delete("pot.dat");
            File.Delete("carrier_density.dat");
            File.Delete("charge_density.dat");
            File.Delete("potential.dat");
            File.Delete("lap.dat");

            Close(converged, max_iterations);
        }

        double dens_diff_lim = 0.1; // the maximum percentage change in the density required for update of V_xc
        double max_vxc_diff = double.MaxValue; // maximum difference for dft potential... if this increases, the dft mixing parameter is reduced
        double min_dens_diff = 0.02; // minimum bound for the required, percentage density difference for updating the dft potential
        double min_vxc_diff = 0.1; // minimum difference in the dft potential for convergence
        double min_alpha = 0.03; // minimum possible value of the dft mixing parameter
        protected override bool Run_Iteration_Routine(IDensity_Solve dens_solv, IPoisson_Solve pois_solv, double tol, int max_iterations)
        {
            dens_solv.Set_DFT_Potential(carrier_density);
            if (!no_dft)
            {
                dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                dens_solv.Set_DFT_Potential(carrier_density);
            }

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
                Band_Data gphi = -1.0 * chem_pot.Laplacian / Physics_Base.q_e - carrier_density.Spin_Summed_Data;
                Set_Edges(gphi);
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, gphi, carrier_density, dens_solv.DFT_diff(carrier_density));
               // chem_pot = pois_solv.Chemical_Potential;

                // Calculate optimal damping parameter, t, (but damped damping....)
                if (t == 0.0)
                    t = t_min;

                t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot / Physics_Base.q_e, x, carrier_density, dopent_density, pois_solv, dens_solv, t_min);
                if (t < 0.0)
                {
                    Console.WriteLine("Iterator has stalled, setting t = 0");
                    t = 0.0;
                }

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
                    if (current_vxc_diff > max_vxc_diff && dens_diff_lim / 2.0 > min_dens_diff || no_dft)
                    {
                        dens_diff_lim /= 2.0;
                        Console.WriteLine("Minimum percentage density difference reduced to " + dens_diff_lim.ToString());
                    }
                    max_vxc_diff = current_vxc_diff;

                    // solution is converged if the density accuracy is better than half the minimum possible value for changing the dft potential
                    if (dens_diff.Max() < min_dens_diff / 2.0 && current_vxc_diff < min_vxc_diff &&  Physics_Base.q_e * x.InfinityNorm() < tol && t != t_min)
                        converged = true;
                }

                // update t for Poisson solver
                pois_solv.T = t;
                chem_pot = chem_pot + t * Physics_Base.q_e * x;

                base.Checkpoint();

                stpwch.Stop();
                Console.WriteLine("Iter = " + count.ToString() + "\tDens = " + dens_diff.Max().ToString("F4") + "\tPot = " + (Physics_Base.q_e * x.InfinityNorm()).ToString("F6") + "\tt = " + t.ToString("F5") + "\ttime = " + stpwch.Elapsed.TotalMinutes.ToString("F"));
                count++;

            //    File.Copy("split_gate.pg6", "split_gate_" + count.ToString("000") + ".pg6");

                // reset the potential if the added potential t * x is too small
                if (converged || count > max_iterations)
                {
                    Console.WriteLine("Maximum potential change at end of iteration was " + (t * Physics_Base.q_e * x.InfinityNorm()).ToString());
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
        void Set_Edges(SpinResolved_Data carrier_density)
        {
            Set_Edges(carrier_density.Spin_Up);
            Set_Edges(carrier_density.Spin_Down);
        }

        void Set_Edges(Band_Data data)
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

            // top and bottom must be set to zero
            for (int i = 0; i < nx_dens; i++)
                for (int j = 0; j < ny_dens; j++)
                {

                    data.vol[0][i, j] = 0.0;
                    data.vol[nz_dens - 1][i, j] = 0.0;
                }
        }
        void Initialise_from_1D(Dictionary<string, object> input_dict)
        {
            // get data from dictionary
            SpinResolved_Data tmp_1d_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
            SpinResolved_Data tmp_1d_dopdens = (SpinResolved_Data)input_dict["Dopent_Density"];
            Band_Data tmp_pot_1d = (Band_Data)input_dict["Chemical_Potential"];

            // calculate where the bottom of the 3D data is
            int offset_min = (int)Math.Round((zmin_dens - zmin_pot) / dz_pot);

            // this is the charge density modulation in the (x, y) plane so, initially, just put it in uniformly
            for (int k = 0; k < nz_dens; k++)
                for (int i = 0; i < nx_dens; i++)
                    for (int j = 0; j < ny_dens; j++)
                    {
                        chem_pot.vol[k][i, j] = tmp_pot_1d.vec[offset_min + j];

                        // do not add anything to the density at the top or bottom of the domain
                        if (k == 0 || k == nz_dens - 1)
                            continue;

                        // carrier data
                        carrier_density.Spin_Up.vol[k][i, j] = tmp_1d_density.Spin_Up.vec[k + offset_min];
                        carrier_density.Spin_Down.vol[k][i, j] = tmp_1d_density.Spin_Down.vec[k + offset_min];
                        
                        // dopent data
                        dopent_density.Spin_Up.vol[k][i, j] = tmp_1d_dopdens.Spin_Up.vec[k + offset_min];
                        dopent_density.Spin_Down.vol[k][i, j] = tmp_1d_dopdens.Spin_Down.vec[k + offset_min];
                    }

            boundary_conditions.Add("surface", (double)input_dict["surface_charge"]);
        }

        protected override void Initialise_from_Hot_Start(Dictionary<string, object> input_dict)
        {
            string spin_up_location = (string)input_dict["spin_up_file"];
            string spin_down_location = (string)input_dict["spin_down_file"];

            // load (spin-resolved) density data
            string[] spin_up_data, spin_down_data;
            try
            {
                spin_up_data = File.ReadAllLines(spin_up_location);
                spin_down_data = File.ReadAllLines(spin_down_location);
            }
            catch (KeyNotFoundException key_e)
            { throw new Exception("Error - Are the file names for the hot start data included in the input file?\n" + key_e.Message); }

            this.carrier_density.Spin_Up = Band_Data.Parse_Band_Data(spin_up_location, spin_up_data, Nx_Dens, Ny_Dens, Nz_Dens);
            this.carrier_density.Spin_Down = Band_Data.Parse_Band_Data(spin_down_location, spin_down_data, Nx_Dens, Ny_Dens, Nz_Dens);

            // and surface charge density
            try { boundary_conditions.Add("surface", double.Parse(File.ReadAllLines((string)input_dict["surface_charge_file"])[0])); }
            catch (KeyNotFoundException key_e) { throw new Exception("Error - Are the file names for the hot start data included in the input file?\n" + key_e.Message); }
        }

    }
}
