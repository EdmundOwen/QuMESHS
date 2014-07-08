using System;
using System.Collections.Generic;
using System.Diagnostics;
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
        double t_damp = 1.0, t_min = 1e-3;
        double edge_min_charge = 1e-5;

        TwoD_ThomasFermiSolver dens_solv;
        TwoD_PoissonSolver_Scaled pois_solv;

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

                        // linearly interpolate with exponential envelope in the transverse direction
                        int j_init = (int)Math.Floor(j * Dz_Dens / Dz_Pot);
                        double y = (i - 0.5 * ny_dens) * dy_dens;
                        //double z = (j - 0.5 * nz_dens) * dz_dens;
                        double envelope = 0.0;//Math.Exp(-100.0 * y * y / ((double)(ny_dens * ny_dens) * dy_dens * dy_dens));// *Math.Exp(-10.0 * z * z / ((double)(nz_dens * nz_dens) * dz_dens * dz_dens)); 
                        carrier_density.Spin_Up.mat[i, j] = envelope * tmp_1d_density.Spin_Up.vec[offset_min + j_init];// +(tmp_1d_density.Spin_Up.vec[offset_min + j_init + 1] - tmp_1d_density.Spin_Up.vec[offset_min + j_init]) * Dz_Dens / Dz_Pot;
                        carrier_density.Spin_Down.mat[i, j] = envelope * tmp_1d_density.Spin_Down.vec[offset_min + j_init];// +(tmp_1d_density.Spin_Down.vec[offset_min + j_init + 1] - tmp_1d_density.Spin_Down.vec[offset_min + j_init]) * Dz_Dens / Dz_Pot;
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

            if (input_dict.ContainsKey("dft")) this.TF_only = !(bool)input_dict["dft"];
            if (input_dict.ContainsKey("TF_only")) this.TF_only = (bool)input_dict["TF_only"];

            // create charge density solver and calculate boundary conditions
            dens_solv = new TwoD_ThomasFermiSolver(this);
            double bottom_V = dens_solv.Get_Chemical_Potential(0.0, zmin_pot, layers) / (Physics_Base.q_e * Physics_Base.energy_V_to_meVpzC);

            // initialise potential solver
            bool with_smoothing = false;
            Get_From_Dictionary<bool>(input_dict, "with_smoothing", ref with_smoothing);
            pois_solv = new TwoD_PoissonSolver_Scaled(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            
            Console.WriteLine("Experimental parameters initialised");
        }

        public override void Run()
        {
            int count = 0;
            bool converged;
            //if (TF_only)
            //{
            t = t_damp;
                converged = false;
                while (!converged)
                {
                    if (count == 0)
                    {
                        // calculate initial potential with the given charge distribution
                        Console.WriteLine("Calculating initial potential grid");
                        pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);
                        chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
                        Console.WriteLine("Initial grid complete");
                        // Get charge rho(phi) (not dopents as these are included as a flexPDE input)
                        dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                    }

                    Stopwatch stpwch = new Stopwatch();
                    stpwch.Start();

                    // Generate the charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi)
                    SpinResolved_Data rho_prime = dens_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                    // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                    //Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, carrier_density.Spin_Summed_Data);
                    Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data, carrier_density.Spin_Summed_Data);
                    chem_pot = pois_solv.Chemical_Potential;

                    // load residual_g_phi and residual_g_x from the output of the Newton step calculation
                    residual_g_phi = Get_FlexPDE_Data("residual_g.dat", "residual_g_phi");
                    residual_g_x = Get_FlexPDE_Data("residual_g.dat", "residual_g_x");

                    // Calculate optimal damping parameter, t, (but damped damping....)
                    t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot, x, carrier_density, dopent_density, pois_solv, dens_solv, t_min);

                    // Check convergence
                    Band_Data g_phi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data;
                    double[] diff = new double[ny_dens * nz_dens];
                    for (int j = 0; j < ny_dens * nz_dens; j++)
                        diff[j] = Math.Abs(g_phi[j]);
                    double convergence = diff.Sum();
                    if (diff.Max() < tol || Math.Max((t * x).mat.Max(), (-t * x).mat.Max()) < 1.0)
                        converged = true;

                    // Recalculate the charge density but for the updated potential rho(phi + t * x)
                    bool edges_fine = false;
                    while (!edges_fine)
                    {
                        edges_fine = true;
                        dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot + t * x);
                        for (int i = 1; i < ny_dens - 1; i++)
                            for (int j = 1; j < nz_dens - 1; j++)
                                //if (i == 1 || j == 1 || i == ny_dens - 2 || j == nz_dens - 2)
                                if (i == 1 || i == ny_dens - 2)
                                    if (Math.Abs(carrier_density.Spin_Summed_Data.mat[i, j]) > edge_min_charge)
                                    {
                                        if (x.mat.Max() > 0)
                                            t = 0.5 * t;
                                        else
                                            t = 2.0 * t;

                                        edges_fine = false;
                                        goto end;
                                    }

                    end:
                        continue;
                    }

                    // update band energy phi_new = phi_old + t * x
                    chem_pot = chem_pot + t * x;
                    pois_solv.T = t;
                    stpwch.Stop();
                    Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString("F") + "\tt = " + t.ToString() + "\ttime = " + stpwch.Elapsed.TotalMinutes.ToString("F"));
                    count++;

                    if ((count - 1) % 10 == 0)
                        File.Copy("split_gate.pg6", "split_gate_TF" + (count - 1).ToString("000") + ".pg6", true);
                    
                }
            //}

                Console.WriteLine("Initial Thomas-Fermi approximation solution has converged!");

            // and then run the DFT solver at the base temperature over a limited range
            TwoD_DFTSolver dft_solv = new TwoD_DFTSolver(this);
            dft_solv.Ymin_Pot = ymin_pot; dft_solv.Dy_Pot = dy_pot;
            dft_solv.Zmin_Pot = zmin_pot; dft_solv.Dz_Pot = dz_pot;
            //TwoD_EffectiveBandSolver dft_solv = new TwoD_EffectiveBandSolver(this);

            // and recalculate the potential for the TF density
            Console.WriteLine("Calculating potential grid for TF density");
            pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);
            chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
            Console.WriteLine("Initial TF potential calculated");
            dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);

            count = 0;
            converged = false;
            while (!converged)
            {
                if (count == 10)
                {
                    Console.WriteLine("ReCalculating with quantum mechanical density");
                    pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);
                    chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
                    Console.WriteLine("Initial QM potential calculated");
                    dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                }

                Stopwatch stpwch = new Stopwatch();
                stpwch.Start();

                // Generate an approximate charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi) using the Thomas-Fermi semi-classical method
                SpinResolved_Data rho_prime = dft_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data, carrier_density.Spin_Summed_Data);
                chem_pot = pois_solv.Chemical_Potential;
                // NOTE: the XC potential must be subtracted as this is not included in the poisson equation
                //Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * pois_solv.Calculate_Laplacian((chem_pot - Physics_Base.Get_XC_Potential(carrier_density)) / Physics_Base.q_e) - carrier_density.Spin_Summed_Data, carrier_density.Spin_Summed_Data);
                //chem_pot = pois_solv.Chemical_Potential + Physics_Base.Get_XC_Potential(carrier_density);

                // load residual_g_phi and residual_g_x from the output of the Newton step calculation
                residual_g_phi = Get_FlexPDE_Data("residual_g.dat", "residual_g_phi");
                residual_g_x = Get_FlexPDE_Data("residual_g.dat", "residual_g_x");

                // Calculate optimal damping parameter, t, (but damped damping....)
                t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot, x, carrier_density, dopent_density, pois_solv, dft_solv, t_min);

                // Check convergence
                Band_Data g_phi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data;
                double[] diff = new double[ny_dens * nz_dens];
                for (int j = 0; j < ny_dens * nz_dens; j++)
                    diff[j] = Math.Abs(g_phi[j]);
                double convergence = diff.Sum();
                if (diff.Max() < tol)
                    converged = true;

                // Recalculate the charge density but for the updated potential rho(phi + t * x)
                bool edges_fine = false;
                while (!edges_fine)
                {
                    edges_fine = true;
                    dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot + t * x);
                    for (int i = 1; i < ny_dens - 1; i++)
                        for (int j = 1; j < nz_dens - 1; j++)
                            //if (i == 1 || j == 1 || i == ny_dens - 2 || j == nz_dens - 2)
                            if (i == 1 || i == ny_dens - 2)
                                if (Math.Abs(carrier_density.Spin_Summed_Data.mat[i, j]) > edge_min_charge)
                                {
                                    t = 2.0 * t;
                                    edges_fine = false;
                                    goto end;
                                }

                    end:
                        continue;
                }

                // update band energy phi_new = phi_old + t * x
                chem_pot = chem_pot + t * x;
                pois_solv.T = t;
                stpwch.Stop();
                Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString("F") + "\tt = " + t.ToString() + "\ttime = " + stpwch.Elapsed.TotalMinutes.ToString("F"));
                count++;

                // reset the potential if the added potential t * x is too small
                if (Math.Max((1.0 * x).mat.Max(), (-1.0 * x).mat.Max()) < 1.0)
                {
                    File.Copy("split_gate.pg6", "split_gate_WF_BRfinal" + (count - 1).ToString("000") + ".pg6", true);
                    Console.WriteLine("Maximum potential change this iteration was " + Math.Max((t * x).mat.Max(), (-t * x).mat.Max()).ToString() + "\nChanging to potential blend method...");
                    break;
                }

                if ((count - 1) % 10 == 0 && count < 3000)
                    File.Copy("split_gate.pg6", "split_gate_WF_BR" + (count - 1).ToString("0000") + ".pg6", true);

                if ((count - 1) % 100 == 0 && count >= 3000)
                    File.Copy("split_gate.pg6", "split_gate_WF_BR" + (count - 1).ToString("0000") + ".pg6", true);
            }

            Console.WriteLine("Bank-Rose iteration complete");
            
            /*
            int no_finaliters = 10;
            for (int iter = 0; iter < no_finaliters; iter++)
            {
                count = 0;
                converged = false;
                while (!converged)
                {
                    Stopwatch stpwch = new Stopwatch();
                    stpwch.Start();

                    // Generate an approximate charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi) using the Thomas-Fermi semi-classical method
                    SpinResolved_Data rho_prime = dft_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                    // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                    Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, carrier_density.Spin_Summed_Data);
                    chem_pot = pois_solv.Chemical_Potential;

                    // load residual_g_phi and residual_g_x from the output of the Newton step calculation
                    residual_g_phi = Get_FlexPDE_Data("residual_g.dat", "residual_g_phi");
                    residual_g_x = Get_FlexPDE_Data("residual_g.dat", "residual_g_x");

                    // Calculate optimal damping parameter, t, (but damped damping....)
                    t = t_damp * Calculate_optimal_t(t / t_damp, chem_pot, x, carrier_density, dopent_density, pois_solv, dft_solv, t_min);

                    // Check convergence
                    Band_Data g_phi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data;
                    double[] diff = new double[ny_dens * nz_dens];
                    for (int j = 0; j < ny_dens * nz_dens; j++)
                        diff[j] = Math.Abs(g_phi[j]);
                    double convergence = diff.Sum();
                    if (diff.Max() < tol)
                        converged = true;

                    // Recalculate the charge density but for the updated potential rho(phi + t * x)
                    bool edges_fine = false;
                    while (!edges_fine)
                    {
                        edges_fine = true;
                        dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot + t * x);
                        for (int i = 1; i < ny_dens - 1; i++)
                            for (int j = 1; j < nz_dens - 1; j++)
                                //if (i == 1 || j == 1 || i == ny_dens - 2 || j == nz_dens - 2)
                                if (i == 1 || i == ny_dens - 2)
                                    if (Math.Abs(carrier_density.Spin_Summed_Data.mat[i, j]) > edge_min_charge)
                                    {
                                        t = 0.5 * t;
                                        edges_fine = false;
                                        goto end;
                                    }

                    end:
                        continue;
                    }

                    // update band energy phi_new = phi_old + t * x
                    chem_pot = chem_pot + t * x;
                    pois_solv.T = t;
                    stpwch.Stop();
                    Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString("F") + "\tt = " + t.ToString() + "\ttime = " + stpwch.Elapsed.TotalMinutes.ToString("F"));
                    count++;

                    // reset the potential if the added potential t * x is too small
                    if (Math.Max((t * x).mat.Max(), (-t * x).mat.Max()) < 0.01)
                    {
                        File.Copy("phi.dat", "phi_process" + iter.ToString("00") + ".dat", true);
                        Console.WriteLine("Maximum potential change this iteration was " + Math.Max((t * x).mat.Max(), (-t * x).mat.Max()).ToString() + "\nResetting potential...");
                        pois_solv.Set_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge);
                        chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);

                        converged = true;
                    }
                }
            }
            
            /*
            t = 0.003;
            count = 0;
            pois_solv.Set_Blend_Boundary_Conditions(top_V, split_V, split_width, bottom_V, surface_charge, t);
            while(!converged)
            {
                Stopwatch stpwch = new Stopwatch();
                stpwch.Start();

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data);
                chem_pot = pois_solv.Chemical_Potential;

                // Check convergence
                Band_Data g_phi = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e) - carrier_density.Spin_Summed_Data;
                double[] diff = new double[ny_dens * nz_dens];
                for (int j = 0; j < ny_dens * nz_dens; j++)
                    diff[j] = Math.Abs(g_phi[j]);
                double convergence = diff.Sum();
                if (diff.Max() < tol)
                    converged = true;

                dft_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);

                // update band energy phi_new = phi_old + t * x
                stpwch.Stop();
                Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString("F") + "\tt = " + t.ToString() + "\ttime = " + stpwch.Elapsed.TotalMinutes.ToString("F"));
                count++;

                if ((count - 1) % 10 == 0)
                    File.Copy("split_gate.pg6", "split_gate_WF_blend" + (count - 1).ToString("000") + ".pg6", true);
            
            }
            */

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

        double Get_FlexPDE_Data(string input_filename, string data_string)
        {
            string[] data = File.ReadAllLines(input_filename);

            for (int i = 0; i < data.Length; i++)
                if (data[i].Contains(data_string))
                    return double.Parse(data[i].Split(' ').Last());

            throw new KeyNotFoundException();
        }

        protected override double calc_vp(double t, Band_Data band_energy, Band_Data x, SpinResolved_Data car_dens, SpinResolved_Data dop_dens, IPoisson_Solve pois_solv, IDensity_Solve dens_solv)
        {
            return base.calc_vp(t, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);
            //return residual_g_phi + t * residual_g_x + base.calc_vp(t, band_energy, x, car_dens, dop_dens, pois_solv, dens_solv);
        }
    }
}
