using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Analysis;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace OneD_ThomasFermiPoisson
{
    public class Experiment : Experiment_Base
    {
        Dictionary<double, double> surface_charge;

        bool fix_bottom_V;
        double t_min = 1e-3;
        double t_damp = 1.0;
        OneD_PoissonSolver pois_solv;

        public new void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            Get_From_Dictionary<double>(input_dict, "dz", ref dz_dens); dz_pot = dz_dens;
            Get_From_Dictionary(input_dict, "nz", ref nz_dens); nz_pot = nz_dens;

            // physics parameters are done by the base method
            base.Initialise(input_dict);

            Initialise_DataClasses(input_dict);

            // check that the size of the domain [(nz_pot-1) * dz_pot] is not larger than the band structure
            if (layers[layers.Length - 1].Zmax - Geom_Tool.Get_Zmin(layers) < (Nz_Pot - 1) * Dz_Pot)
                throw new Exception("Error - the band structure provided is smaller than the simulation domain!\nUse nz = " + (int)Math.Ceiling((layers[layers.Length - 1].Zmax - Geom_Tool.Get_Zmin(layers)) / Dz_Pot) + " instead");
            // and check that the top of the domain is the surface (unless this check is overloaded)
            bool surface_override = false; Get_From_Dictionary<bool>(input_dict, "surface_check", ref surface_override, true);
            if (!surface_override && Geom_Tool.Find_Layer_Below_Surface(layers).Zmax - Geom_Tool.Get_Zmin(layers) - Nz_Pot * Dz_Pot != 0.0)
                throw new Exception("Error - the top of the domain is not the surface!\nUse the input \"surface_check\" to override");

            // calculate the top and bottom of the domain
            input_dict["zmin"] = Geom_Tool.Get_Zmin(layers);
            input_dict["zmax"] = Geom_Tool.Get_Zmin(layers) + dz_pot * nz_pot;

            // and split gate dimensions
            device_dimensions.Add("bottom_position", (double)input_dict["zmin"]);
            device_dimensions.Add("top_position", (double)input_dict["zmax"]);

            // check whether the bottom should be fixed
            if (input_dict.ContainsKey("bottom_V"))
                fix_bottom_V = true;

            // initialise the dictionary which contains the surface charges at each temperature
            surface_charge = new Dictionary<double, double>();

            // Initialise potential solver
            pois_solv = new OneD_PoissonSolver(this, using_flexPDE, input_dict);
            
            Console.WriteLine("Experimental parameters initialised");
        }

        protected override void Initialise_DataClasses(Dictionary<string, object> input_dict)
        {
            // try to get and the carrier and dopent density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.carrier_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
            else this.carrier_density = new SpinResolved_Data(nz_pot);
            if (input_dict.ContainsKey("SpinResolved_Dopent_Density")) this.dopent_density = (SpinResolved_Data)input_dict["SpinResolved_Dopent_Density"];
            else this.dopent_density = new SpinResolved_Data(nz_pot);
            // and instantiate their derivatives
            carrier_density_deriv = new SpinResolved_Data(nz_pot);
            dopent_density_deriv = new SpinResolved_Data(nz_pot);

            // and finally, try to get the chemical potential from the dictionary...
            if (input_dict.ContainsKey("Chemical_Potential")) this.chem_pot = new Band_Data((DoubleVector)input_dict["Chemical_Potential"]); else chem_pot = new Band_Data(nz_pot, 0.0);
        }

        public override void Run()
        {
            // get temperatures to run the experiment at
            double[] run_temps = Freeze_Out_Temperatures();
            
            // run experiment using Thomas-Fermi solver
            for (int i = 0; i < run_temps.Length; i++)
            {
                double current_temperature = run_temps[i];
                this.temperature = current_temperature;

                OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(this, Dz_Pot, Zmin_Pot, Nz_Pot);
                dens_solv.DFT_Mixing_Parameter = 0.0;

                if (!Geom_Tool.GetLayer(layers, zmin_pot).Dopents_Frozen_Out(current_temperature) && !fix_bottom_V)
                    boundary_conditions["bottom_V"] = dens_solv.Get_Chemical_Potential(zmin_pot, layers) / (Physics_Base.q_e * Physics_Base.energy_V_to_meVpzC);

                Run_Iteration_Routine(dens_solv, tol);

                // save the surface charge for this temperature
                surface_charge.Add(current_temperature, pois_solv.Get_Surface_Charge(chem_pot, layers));
                pois_solv.Reset();
            }

            if (!TF_only)
            {
                // and then run the DFT solver at the base temperature
                OneD_DFTSolver dft_solv = new OneD_DFTSolver(this);
                dft_solv.DFT_Mixing_Parameter = 0.0;                 //NOTE: This method doesn't mix in the DFT potential in this way (DFT is still implemented)
                dft_solv.Zmin_Pot = zmin_pot; dft_solv.Dz_Pot = dz_pot;
                Console.WriteLine("Starting DFT calculation");
                Run_Iteration_Routine(dft_solv, tol);

                pois_solv.Reset();

                // initialise output solvers
                OneD_ThomasFermiSolver final_dens_solv = new OneD_ThomasFermiSolver(this, Dz_Pot, Zmin_Pot, Nz_Pot);

                (Input_Band_Structure.Get_BandStructure_Grid(layers, dz_pot, nz_pot, zmin_pot) - chem_pot).Save_Data("potential" + output_suffix);
            }

            double tot_dens = (from val in carrier_density.Spin_Summed_Data.vec
                               where val < 0.0
                               select -1.0e14 * val * dz_dens / Physics_Base.q_e).ToArray().Sum();
            Console.WriteLine("Carrier density at heterostructure interface: \t" + tot_dens.ToString("e3") + " cm^-2");

            // there is no iteration timeout for the 1D solver so if it gets to this point the solution will definitely have converged
            Close(true, int.MaxValue);
        }

        bool Run_Iteration_Routine(IDensity_Solve dens_solv, double pot_lim)
        {
            return Run_Iteration_Routine(dens_solv, pot_lim, int.MaxValue);
        }

        bool Run_Iteration_Routine(IDensity_Solve dens_solv, double pot_lim, int max_count)
        {
            // calculate initial potential with the given charge distribution
            Console.WriteLine("Calculating initial potential grid");
            pois_solv.Initiate_Poisson_Solver(device_dimensions, boundary_conditions);
            chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data);
            Console.WriteLine("Initial grid complete");
            dens_solv.Set_DFT_Potential(carrier_density);
            dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
            dens_solv.Set_DFT_Potential(carrier_density);

            int count = 0;
            t = 1.0;
            bool converged = false;
            while (!converged)
            {
                Band_Data dens_old = carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data;

                // Get charge rho(phi)
                dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);

                // Generate charge-dependent part of the Jacobian, g'(phi) = -d(eps * d( )) - rho'(phi)
                SpinResolved_Data rho_prime = dens_solv.Get_ChargeDensity_Deriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                // Calculate Laplacian operating on the given band energy, d(eps * d(phi))
                Band_Data g_phi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data - dopent_density.Spin_Summed_Data;
                g_phi[0] = 0.0; g_phi[g_phi.Length - 1] = 0.0;
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, g_phi);

                // Calculate optimal damping parameter, t
                t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot, x, carrier_density, dopent_density, pois_solv, dens_solv, t_min);

                // Check convergence
                Band_Data dens_spin_summed = carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data;
                Band_Data dens_diff = dens_spin_summed - dens_old;
                double dens_min = Math.Max(Math.Abs(dens_spin_summed.Max()), Math.Abs(dens_spin_summed.Min()));
                for (int i = 0; i < dens_diff.Length; i++)
                    // only calculate density difference for densities more than 1% of the maximum value
                    if (Math.Abs(dens_spin_summed[i]) > 0.01 * dens_min)
                        dens_diff[i] = Math.Abs(dens_diff[i] / dens_spin_summed[i]);
                    else
                        dens_diff[i] = 0.0;

                double[] diff = new double[Nz_Pot];
                for (int j = 0; j < nz_pot; j++)
                    diff[j] = Math.Abs(g_phi.vec[j]);
                double convergence = diff.Sum();
                if (diff.Max() < tol)
                    converged = true;

                // update band energy phi_new = phi_old + t * x
               Console.WriteLine("Iter = " + count.ToString() + "\tDens conv = " + dens_diff.Max().ToString("F4") + "\tt = " + t.ToString());
                chem_pot = chem_pot + t * x;
                count++;

                // reset the potential if the added potential t * x is too small
                if (converged || count > max_count)
                {
                    Console.WriteLine("Maximum potential change at end of iteration was " + Math.Max(t * x.Max(), (-t * x).Max()).ToString());
                    break;
                }
            }

            Console.WriteLine("Iteration complete");

            return converged;
        }

        private void Run_Test()
        {
            throw new NotImplementedException();
            /*
            StreamWriter sw_dens = new StreamWriter("test_dens.dat");
            StreamWriter sw_e = new StreamWriter("test_ener.dat");

            Console.WriteLine("top_gate\t\tN_top\t\tN_bottom");
            int top_index = (int)Math.Floor((- zmin_pot + zmin_dens) / dz_pot + nz_dens / 2);
            int bottom_index = (int)Math.Floor((- zmin_pot + zmin_dens) / dz_pot);
            int test_no = 50;
            for (int i = 0; i < test_no; i++)
            {
                top_V = -0.02 * (double)i;

                // and then run the DFT solver at the base temperature
                OneD_DFTSolver dft_solv = new OneD_DFTSolver(this);
                dft_solv.DFT_Mixing_Parameter = 0.0;                 //NOTE: This method doesn't mix in the DFT potential in this way (DFT is still implemented)
                dft_solv.Zmin_Pot = zmin_pot; dft_solv.Dz_Pot = dz_pot;
                Console.WriteLine("Starting DFT calculation");
                Run_Iteration_Routine(dft_solv, tol);

                double top = 0.0;
                double bottom = 0.0;
                Band_Data car_dens_spin_summed = carrier_density.Spin_Summed_Data;
                for (int j= 0; j < nz_dens / 2; j++)
                {
                    top += car_dens_spin_summed.vec[top_index + j];
                    bottom += car_dens_spin_summed.vec[bottom_index + j];
                }
                top *= -1.0e14 * dz_dens / Physics_Base.q_e;
                bottom *= -1.0e14 * dz_dens / Physics_Base.q_e;

                Console.WriteLine(top_V.ToString() + "\t\t" + top.ToString("e3") + "\t\t" + bottom.ToString("e3"));
                sw_dens.WriteLine(top.ToString("e3") + "\t" + bottom.ToString("e3"));
                DoubleVector energies = dft_solv.Get_EnergyLevels(layers, ref carrier_density, chem_pot);
                sw_e.WriteLine(energies[0].ToString() + "\t" + energies[1].ToString() + "\t" + energies[2].ToString() + "\t" + energies[3].ToString());

                File.Copy("tmp_pot.dat", "pot_" + i.ToString("00"), true);
                File.Copy("tmp_dens.dat", "dens_" + i.ToString("00"), true);
            }

            sw_dens.Close();
            sw_e.Close();*/
        }

        /// <summary>
        /// returns the precalculated surface charge for this temperature
        /// throws an error if this isn't exactly the temperature a simulation was calculated at
        /// </summary>
        public double Surface_Charge(double temperature)
        {
            if (!surface_charge.ContainsKey(temperature))
                throw new KeyNotFoundException("Error - the charge density was not calculated for T = " + temperature.ToString() + "K");

            return surface_charge[temperature];
        }

        protected override void Initialise_from_Hot_Start(Dictionary<string, object> input_dict)
        {
            throw new NotImplementedException();
        }
    }
}