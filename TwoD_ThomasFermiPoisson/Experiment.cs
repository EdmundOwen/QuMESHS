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
        double residual_g_phi, residual_g_x;

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

            // and split gate dimensions
            Get_From_Dictionary<double>(input_dict, "split_width", ref split_width);

            // try to get the density from the dictionary... they probably won't be there and if not... make them
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
                this.dopent_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));

                SpinResolved_Data tmp_1d_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
                SpinResolved_Data tmp_1d_dopdens = (SpinResolved_Data)input_dict["Dopent_Density"];

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
                        dopent_density.Spin_Up.mat[i, j] = tmp_1d_dopdens.Spin_Up.vec[offset_min + j_init] + (tmp_1d_dopdens.Spin_Up.vec[offset_min + j_init + 1] - tmp_1d_dopdens.Spin_Up.vec[offset_min + j_init]) * Dz_Dens / Dz_Pot;
                        dopent_density.Spin_Down.mat[i, j] = tmp_1d_dopdens.Spin_Down.vec[offset_min + j_init] + (tmp_1d_dopdens.Spin_Down.vec[offset_min + j_init + 1] - tmp_1d_dopdens.Spin_Down.vec[offset_min + j_init]) * Dz_Dens / Dz_Pot;
                    }

                Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            }
            else
            {
                this.carrier_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));
                this.dopent_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));
                Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            }

            // and instantiate their derivatives
            carrier_density_deriv = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));
            dopent_density_deriv = new SpinResolved_Data(new Band_Data(new DoubleMatrix(ny_dens, nz_dens)), new Band_Data(new DoubleMatrix(ny_dens, nz_dens)));
            
            // and do the same for the chemical potential
            if (input_dict.ContainsKey("Chemical_Potential"))
            {
                chem_pot = new Band_Data(new DoubleMatrix(ny_dens, nz_dens));
                Band_Data tmp_pot_1d = (Band_Data)input_dict["Chemical_Potential"];

                int offset_min = (int)Math.Round((zmin_dens - zmin_pot) / dz_pot);
                for (int i = 0; i < ny_dens; i++)
                    for (int j = 0; j < nz_dens; j++)
                    {
                        chem_pot.mat[i, j] = tmp_pot_1d.vec[offset_min + j - 1];
                    }
            }
            else
                chem_pot = new Band_Data(new DoubleMatrix(ny_dens, nz_dens));

            if (input_dict.ContainsKey("dft")) this.TF_only = !bool.Parse((string)input_dict["dft"]);
            if (input_dict.ContainsKey("TF_only")) this.TF_only = bool.Parse((string)input_dict["TF_only"]);

            // create charge density solver and calculate boundary conditions
            dens_solv = new TwoD_ThomasFermiSolver(this);
            double bottom_V = dens_solv.Get_Chemical_Potential(0.0, zmin_pot, layers) / (Physics_Base.q_e * Physics_Base.energy_V_to_meVpzC);

            // initialise potential solver
            pois_solv = new TwoD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);

            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            int count = 0;
            bool converged;
            
            if (TF_only)
            {
                t = 1.0;
                converged = false;
                while (!converged)
                {
                    // Get charge rho(phi) (not dopents as these are included as a flexPDE input)
                    dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                    Band_Data charge_dens_old = carrier_density.Spin_Summed_Data;

                    // Calculate Laplacian operating on the given band energy d(eps * d(phi))
                    Band_Data tmp_g = pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e);

                    // Generate Jacobian g'(phi) = d(eps * d( )) + rho'(phi)
                    SpinResolved_Data rho_prime = dens_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                    // Solve stepping equation to find raw Newton iteration step, g'(phi) x = g(phi)
                    Band_Data g_u = tmp_g + charge_dens_old;
                    Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * g_u);

                    // Calculate optimal damping parameter, t
                    t = Calculate_optimal_t(t, chem_pot, x, carrier_density, dopent_density, pois_solv, dens_solv);

                    // Check convergence
                    double[] diff = new double[ny_pot * nz_pot];
                    for (int j = 0; j < ny_pot * nz_pot; j++)
                        diff[j] = Math.Abs(g_u.vec[j]);
                    double convergence = diff.Sum();
                    if (diff.Max() < tol)
                        converged = true;

                    // update band energy phi_new = phi_old + t * x
                    Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString() + "\tt = " + t.ToString());
                    chem_pot = chem_pot + t * x;
                    count++;
                }

                pois_solv.Reset();
            }

            // and then run the DFT solver at the base temperature over a limited range
            TwoD_DFTSolver dft_solv = new TwoD_DFTSolver(this);
            dft_solv.Ymin_Pot = ymin_pot; dft_solv.Dy_Pot = dy_pot;
            dft_solv.Zmin_Pot = zmin_pot; dft_solv.Dz_Pot = dz_pot;
            //TwoD_EffectiveBandSolver dft_solv = new TwoD_EffectiveBandSolver(this);

            // calculate initial potential with the given charge distribution
            chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);

            count = 0;
            t = 1.0;
            converged = false;
            while (!converged)
            {
                // Get charge rho(phi) (not dopents as these are included as a flexPDE input)
                dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                Band_Data charge_dens_old = carrier_density.Spin_Summed_Data;

                // Calculate Laplacian operating on the given band energy d(eps * d(phi))
                Band_Data tmp_g = pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e);

                // Generate Jacobian g'(phi) = d(eps * d( )) + rho'(phi)
                SpinResolved_Data rho_prime = dft_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                Band_Data g_u = tmp_g + charge_dens_old;
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * g_u);

                // load residual_g_phi and residual_g_x from the output of the Newton step calculation

                // Calculate optimal damping parameter, t
                t = Calculate_optimal_t(t, chem_pot, x, carrier_density, dopent_density, pois_solv, dft_solv);
                pois_solv.Update_Potential(t);

                // Check convergence
                double[] diff = new double[ny_dens * nz_dens];
                for (int j = 0; j < ny_dens * nz_dens; j++)
                    diff[j] = Math.Abs(g_u[j]);
                double convergence = diff.Sum();
                if (diff.Max() < tol)
                    converged = true;

                // update band energy phi_new = phi_old + t * x
                Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString() + "\tt = " + t.ToString());
                chem_pot = chem_pot + t * x;
                count++;
            }

            // initialise output solvers
            TwoD_ThomasFermiSolver final_dens_solv = new TwoD_ThomasFermiSolver(this);
            TwoD_PoissonSolver final_pois_solv = new TwoD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            carrier_density.Spin_Summed_Data.Save_Data("dens_2D_raw.dat");
            carrier_density.Spin_Up.Save_Data("dens_2D_up_raw.dat");
            carrier_density.Spin_Down.Save_Data("dens_2D_down_raw.dat");

            // save surface charge
            StreamWriter sw = new StreamWriter("surface_charge.dat"); sw.WriteLine(surface_charge.ToString()); sw.Close();
            // save eigen-energies
            DoubleVector energies = dft_solv.Get_EnergyLevels(layers, carrier_density, chem_pot);
            StreamWriter sw_e = new StreamWriter("energies.dat");
            for (int i = 0; i < energies.Length; i++)
                sw_e.WriteLine(energies[i]);
            sw_e.Close();

            final_dens_solv.Output(carrier_density, "carrier_density.dat");
            final_dens_solv.Output(carrier_density - dft_solv.Get_ChargeDensity(layers, carrier_density, dopent_density, chem_pot), "density_error.dat");
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dy_dens, dz_dens, ny_dens, nz_dens, ymin_dens, zmin_dens) - chem_pot, "potential.dat");
        }

        void Tmp_Print(Band_Data data, string file)
        {
            StreamWriter sw = new StreamWriter(file);
            for (int i = 0; i < ny_dens; i++)
                for (int j = 0; j < nz_dens; j++)
                {
                    sw.Write(data.mat[i, j].ToString() + '\t');
                    if (j == nz_dens - 1)
                        sw.WriteLine();
                }
            sw.Close();
        }

        protected override double calc_vp(double t, Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            return residual_g_phi + t * residual_g_x + base.calc_vp(t, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);
        }
    }
}
