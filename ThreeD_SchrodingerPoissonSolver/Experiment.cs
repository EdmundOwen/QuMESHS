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
        double well_depth;
        double t_damp = 1.0, t_min = 1e-3;

        DoubleVector dens_1d;
        TwoD_ThomasFermiSolver dens_solv;
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
            Get_From_Dictionary<double>(input_dict, "surface_charge", ref surface_charge);
            Get_From_Dictionary<double>(input_dict, "bottom_bc", ref bottom_V);

            // and split gate dimensions
            Get_From_Dictionary<double>(input_dict, "split_width", ref split_width);
            Get_From_Dictionary<double>(input_dict, "split_length", ref split_length);
            Get_From_Dictionary<double>(input_dict, "top_length", ref top_length);

            // and well depth
            well_depth = layers[1].Zmax - 1;

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density"))
            {
                this.carrier_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(nx_dens, ny_dens)), new Band_Data(new DoubleMatrix(nx_dens, ny_dens)));
                SpinResolved_Data tmp_charge_1d_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
                dens_1d = new DoubleVector(nz_dens);
                int z_offset = (int)Math.Abs((Zmin_Pot - Zmin_Dens) / Dz_Pot);

                for (int i = 0; i < nz_dens; i++)
                    dens_1d[i] = tmp_charge_1d_density.Spin_Summed_Data.vec[z_offset + i];

                // this is the charge density modulation in the (x, y) plane so, initially, it is set to one
                for (int i = 0; i < nx_dens; i++)
                    for (int j = 0; j < ny_dens; j++)
                    {
                        carrier_density.Spin_Up.mat[i, j] = 0.5;
                        carrier_density.Spin_Down.mat[i, j] = 0.5;
                    }
            }
            else
                this.carrier_density = new SpinResolved_Data(new Band_Data(new DoubleMatrix(nx_dens, ny_dens)), new Band_Data(new DoubleMatrix(nx_dens, ny_dens)));

            if (input_dict.ContainsKey("dft")) this.TF_only = !bool.Parse((string)input_dict["dft"]);
            if (input_dict.ContainsKey("TF_only")) this.TF_only = bool.Parse((string)input_dict["TF_only"]);

            // create charge density solver and calculate boundary conditions
            dens_solv = new TwoD_ThomasFermiSolver(this);

            // initialise potential solver
            pois_solv = new ThreeD_PoissonSolver(this, dens_1d, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            pois_solv.Set_Boundary_Conditions(top_V, split_V, top_length, split_width, split_length, bottom_V, surface_charge);

            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            if (!hot_start)
                Run_Iteration_Routine(dens_solv, 1.0);

            bool converged = false;
            int no_runs = 10;
            while (!converged)
            {
                converged = Run_Iteration_Routine(dens_solv, tol, no_runs);
                no_runs += 10;
            }

            // initialise output solvers
           // TwoD_ThomasFermiSolver final_dens_solv = new TwoD_ThomasFermiSolver(temperature, dx_dens, dy_dens, nx_dens, ny_dens, xmin_dens, ymin_dens);
            ThreeD_PoissonSolver final_pois_solv = new ThreeD_PoissonSolver(this, dens_1d, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            carrier_density.Spin_Summed_Data.Save_Data("dens_2D.dat");
            carrier_density.Spin_Up.Save_Data("dens_2D_up.dat");
            carrier_density.Spin_Down.Save_Data("dens_2D_down.dat");

        //    final_dens_solv.Output(carrier_density, "carrier_density.dat");
       //     final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, nx_dens, ny_dens, xmin_dens, ymin_dens) - chem_pot, "potential.dat");

            throw new NotImplementedException();
        }

        bool Run_Iteration_Routine(IDensity_Solve dens_solv, double pot_lim)
        {
            return Run_Iteration_Routine(dens_solv, pot_lim, int.MaxValue);
        }

        bool Run_Iteration_Routine(IDensity_Solve dens_solv, double pot_lim, int max_count)
        {
            // calculate initial potential with the given charge distribution
            Console.WriteLine("Calculating initial potential grid");
            pois_solv.Set_Boundary_Conditions(top_V, split_V, top_length, split_width, split_length, bottom_V, surface_charge);
            chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
            Console.WriteLine("Initial grid complete");
            // Get charge rho(phi) (not dopents as these are included as a flexPDE input)
            dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);

            int count = 0;
            bool converged = false;
            while (!converged)
            {
                Stopwatch stpwch = new Stopwatch();
                stpwch.Start();

                // Generate an approximate charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi) using the Thomas-Fermi semi-classical method
                SpinResolved_Data rho_prime = dens_solv.Get_ChargeDensity_Deriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data, carrier_density.Spin_Summed_Data);
                chem_pot = pois_solv.Chemical_Potential;

                // Calculate optimal damping parameter, t, (but damped damping....)
                t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot, x, carrier_density, dopent_density, pois_solv, dens_solv, t_min);
                if (t < 0.0)
                {
                    Console.WriteLine("Iterator has stalled, resetting potential");
                    pois_solv.Set_Boundary_Conditions(top_V, split_V, top_length, split_width, split_length, bottom_V, surface_charge);
                    chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
                    dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                    Console.WriteLine("Potential reset!");
                    t = t_min;
                }

                // Check convergence
                Band_Data g_phi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data;
                double[] diff = new double[ny_dens * nz_dens];
                for (int j = 0; j < ny_dens * nz_dens; j++)
                    diff[j] = Math.Abs(g_phi[j]);
                double convergence = diff.Sum();
                if (Math.Max(t * x.mat.Max(), (-t * x).mat.Max()) < pot_lim)
                    converged = true;

                // Recalculate the charge density but for the updated potential rho(phi + t * x)
                dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot + t * x);

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
                    Console.WriteLine("Maximum potential change at end of iteration was " + Math.Max(t * x.mat.Max(), (-t * x).mat.Max()).ToString());
                    break;
                }
            }

            Console.WriteLine("Iteration complete");

            return converged;
        }

    }
}
