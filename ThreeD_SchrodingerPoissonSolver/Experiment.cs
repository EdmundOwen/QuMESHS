using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using CenterSpace.NMath.Core;
using TwoD_ThomasFermiPoisson;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace ThreeD_SchrodingerPoissonSolver
{
    public class Experiment : Experiment_Base
    {
        double top_V, split_V, surface_charge;
        double split_width, split_length;
        double well_depth;

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
            dens_solv = new TwoD_ThomasFermiSolver(temperature, dx_dens, dy_dens, nx_dens, ny_dens, xmin_dens, ymin_dens);

            // initialise potential solver
            pois_solv = new ThreeD_PoissonSolver(this, dens_1d, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, split_length, bottom_V, surface_charge);

            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            int count = 0;/*
            while (!dens_solv.Converged)
            {
                Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + temperature.ToString() + "\tConvergence factor: " + dens_solv.Convergence_Factor.ToString());

                // solve the chemical potential for the given charge density
                chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);

                // find the density for this new chemical potential and blend
                SpinResolved_Data new_density = dens_solv.Get_ChargeDensity(layers, carrier_density, chem_pot, well_depth);
                dens_solv.Blend(ref carrier_density, new_density, alpha, tol);

                //pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                count++;
            }*/

            if (TF_only)
                return;

            // and then run the DFT solver at the base temperature over a limited range and with a reduced mixing constant
            TwoD_DFTSolver dft_solv = new TwoD_DFTSolver(this.Temperature, 1.0, Dx_Dens, Dy_Dens, 1, Nx_Dens, Ny_Dens, double.MaxValue, Xmin_Dens, Ymin_Dens);

            pois_solv.Reset();

            count = 0;
            t = 1.0;
            bool converged = false;
            while (!converged)
            {
                // Get charge rho(phi)
                dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                Band_Data charge_dens_old = carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data;

                // Calculate Laplacian operating on the given band energy d(eps * d(phi))
                Band_Data tmp_g = pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e);

                // Generate Jacobian g'(phi) = d(eps * d( )) + rho'(phi)
                dens_solv.Get_ChargeDensityDeriv(layers, ref carrier_density_deriv, ref dopent_density_deriv, chem_pot);
                Band_Data rho_prime = carrier_density_deriv.Spin_Summed_Data + dopent_density_deriv.Spin_Summed_Data;

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = g(phi)
                Band_Data g_u = tmp_g + charge_dens_old;
                Band_Data x = Physics_Base.q_e * pois_solv.Calculate_Newton_Step(rho_prime, -g_u.vec);

                // Calculate optimal damping parameter, t
                t = Calculate_optimal_t(t, chem_pot, x, carrier_density, dopent_density, pois_solv, dft_solv);

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

            // initialise output solvers
            TwoD_ThomasFermiSolver final_dens_solv = new TwoD_ThomasFermiSolver(temperature, dx_dens, dy_dens, nx_dens, ny_dens, xmin_dens, ymin_dens);
            ThreeD_PoissonSolver final_pois_solv = new ThreeD_PoissonSolver(this, dens_1d, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            carrier_density.Spin_Summed_Data.Save_Data("dens_2D.dat");
            carrier_density.Spin_Up.Save_Data("dens_2D_up.dat");
            carrier_density.Spin_Down.Save_Data("dens_2D_down.dat");

            final_dens_solv.Output(carrier_density, "carrier_density.dat");
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dx_dens, dy_dens, nx_dens, ny_dens, xmin_dens, ymin_dens) - chem_pot, "potential.dat");

            throw new NotImplementedException();
        }
    }
}
