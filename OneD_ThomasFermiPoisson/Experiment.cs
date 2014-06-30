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
        double top_V = 0.0;
        double top_V_cooldown;
        double t_min = 1e-6;

        public new void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            Get_From_Dictionary<double>(input_dict, "dz", ref dz_dens); dz_pot = dz_dens;
            Get_From_Dictionary(input_dict, "nz", ref nz_dens); nz_pot = nz_dens;

            // but try to get the specific values
            Get_From_Dictionary<double>(input_dict, "dz_dens", ref dz_dens, true);
            Get_From_Dictionary<double>(input_dict, "dz_pot", ref dz_pot, true);
            Get_From_Dictionary(input_dict, "nz_dens", ref nz_dens, true);
            Get_From_Dictionary(input_dict, "nz_pot", ref nz_pot, true);

            // physics parameters are done by the base method
            base.Initialise(input_dict);

            // check that the size of the domain [(nz_pot-1) * dz_pot] is not larger than the band structure
            if (layers[layers.Length - 1].Zmax - layers[1].Zmin < (Nz_Pot - 1) * Dz_Pot)
                throw new Exception("Error - the band structure provided is smaller than the simulation domain!\nUse nz = " + (int)Math.Ceiling((layers[layers.Length - 1].Zmax - layers[1].Zmin) / Dz_Pot) + " instead");

            // see whether there's a top gate voltage for determining top_bc
            Get_From_Dictionary<double>(input_dict, "top_V", ref top_V, true);
            // set the cooldown voltage as this and if there is an extra parameter to allow for biased cooldowns then use this
            top_V_cooldown = top_V;
            Get_From_Dictionary<double>(input_dict, "top_V_cooldown", ref top_V_cooldown, true);

            // try to get and the carrier and dopent density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.carrier_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
            else this.carrier_density = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_pot)), new Band_Data(new DoubleVector(nz_pot)));
            if (input_dict.ContainsKey("SpinResolved_Dopent_Density")) this.dopent_density = (SpinResolved_Data)input_dict["SpinResolved_Dopent_Density"];
            else this.dopent_density = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_pot)), new Band_Data(new DoubleVector(nz_pot)));
            // and instantiate their derivatives
            carrier_density_deriv = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_pot)), new Band_Data(new DoubleVector(nz_pot)));
            dopent_density_deriv = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_pot)), new Band_Data(new DoubleVector(nz_pot)));

            // and finally, try to get the chemical potential from the dictionary...
            if (input_dict.ContainsKey("Chemical_Potential")) this.chem_pot = new Band_Data((DoubleVector)input_dict["Chemical_Potential"]); else chem_pot = new Band_Data(new DoubleVector(nz_pot));

            if (input_dict.ContainsKey("dft")) this.TF_only = !(bool)input_dict["dft"];
            if (input_dict.ContainsKey("TF_only")) this.TF_only = (bool)input_dict["TF_only"];
        }

        int count = 0;
        public override void Run()
        {
            // get temperatures to run the experiment at
            double[] run_temps = Freeze_Out_Temperatures();

            // initialise potential solver
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            // set the boundary conditions
            pois_solv.Set_Boundary_Conditions(layers, top_V_cooldown, bottom_V, Geom_Tool.Get_Zmin(layers) + dz_pot * nz_pot, Geom_Tool.Get_Zmin(layers));
            // initialise the chemical potential as the solution with zero density
            chem_pot = pois_solv.Get_Chemical_Potential(carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data);

            for (int i = 0; i < run_temps.Length; i++)
            {
                double current_temperature = run_temps[i];

                OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(current_temperature, dz_pot, nz_pot, zmin_pot);
                if (!Geom_Tool.GetLayer(layers, zmin_pot).Dopents_Frozen_Out(current_temperature))
                    this.bottom_V = dens_solv.Get_Chemical_Potential(zmin_pot, layers) / (Physics_Base.q_e * Physics_Base.energy_V_to_meVpzC);
                // reset the boundary conditions
                pois_solv.Set_Boundary_Conditions(layers, top_V_cooldown, bottom_V, Geom_Tool.Get_Zmin(layers) + dz_pot * nz_pot, Geom_Tool.Get_Zmin(layers));
                
                count = 0;
                t = 1.0;
                bool converged = false;
                while (!converged)
                {
                    // Get charge rho(phi)
                    dens_solv.Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
                    Band_Data charge_dens_old = carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data;

                    // Calculate Laplacian operating on the given band energy, - d(eps * d(phi))
                    Band_Data tmp_g = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e);

                    // Generate charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi)
                    SpinResolved_Data rho_prime = -1.0 * dens_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                    // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                    Band_Data g_phi = tmp_g - charge_dens_old;
                    Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * g_phi);

                    // Calculate optimal damping parameter, t
                    t = Calculate_optimal_t(t, chem_pot, x, carrier_density, dopent_density, pois_solv, dens_solv, t_min);

                    // Check convergence
                    double[] diff = new double[Nz_Pot];
                    for (int j = 0; j < nz_pot; j++)
                        diff[j] = Math.Abs(g_phi.vec[j]);
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

            if (TF_only)
            {
                double tot_dens_TF = (from val in carrier_density.Spin_Summed_Data.vec
                                   where val < 0.0
                                   select -1.0e14 * val * dz_dens / Physics_Base.q_e).ToArray().Sum();
                Console.WriteLine("Carrier density at heterostructure interface: \t" + tot_dens_TF.ToString("e3") + " cm^-2");
                 return;
            }

            // and then run the DFT solver at the base temperature over a limited range
            OneD_DFTSolver dft_solv = new OneD_DFTSolver(temperature, dz_dens, nz_dens, zmin_dens);
            dft_solv.Zmin_Pot = zmin_pot; dft_solv.Dz_Pot = dz_pot;
            Console.WriteLine("Starting DFT calculation");

            // reset boundary conditions
            pois_solv.Set_Boundary_Conditions(layers, top_V, bottom_V, Geom_Tool.Get_Zmin(layers) + dz_pot * nz_pot, Geom_Tool.Get_Zmin(layers));

            count = 0;
            t = 1.0;
            bool dft_converged = false;
            while (!dft_converged)
            {
                // Get charge rho(phi)
                dft_solv.Get_ChargeDensity(layers, ref carrier_density, chem_pot);
                Band_Data charge_dens_old = carrier_density.Spin_Summed_Data + dopent_density.Spin_Summed_Data;
                
                // Calculate Laplacian operating on the given band energy, - d(eps * d(phi))
                Band_Data tmp_g = -1.0 * pois_solv.Calculate_Laplacian(chem_pot / Physics_Base.q_e);

                // Generate an approximate charge-dependent part of the Jacobian, g'(phi) = - d(eps * d( )) - rho'(phi) using the Thomas-Fermi semi-classical method
                SpinResolved_Data rho_prime = -1.0 * dft_solv.Get_ChargeDensityDeriv(layers, carrier_density_deriv, dopent_density_deriv, chem_pot);

                // Solve stepping equation to find raw Newton iteration step, g'(phi) x = - g(phi)
                Band_Data g_phi = tmp_g - charge_dens_old;
                Band_Data x = pois_solv.Calculate_Newton_Step(rho_prime, -1.0 * g_phi);

                // Calculate optimal damping parameter, t
                t = Calculate_optimal_t(t, chem_pot, x, carrier_density, dopent_density, pois_solv, dft_solv, t_min);

                // Check convergence
                double[] diff = new double[Nz_Pot];
                for (int j = 0; j < nz_pot; j++)
                    diff[j] = Math.Abs(g_phi.vec[j]);
                double convergence = diff.Sum();
                if (diff.Max() < tol)
                    dft_converged = true;

                // update band energy phi_new = phi_old + t * x
                Console.WriteLine("Iter = " + count.ToString() + "\tConv = " + convergence.ToString() + "\tt = " + t.ToString());
                chem_pot = chem_pot + t * x;
                count++;
            }

            pois_solv.Reset();

            // initialise output solvers
            OneD_ThomasFermiSolver final_dens_solv = new OneD_ThomasFermiSolver(temperature, dz_dens, nz_dens, zmin_dens);
            OneD_PoissonSolver final_pois_solv = new OneD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            final_dens_solv.Output(carrier_density, "carrier_density.dat", false);

            carrier_density.Spin_Summed_Data.Save_1D_Data("dens_1D.dat", dz_dens, zmin_dens);
            carrier_density.Spin_Up.Save_1D_Data("dens_1D_up.dat", dz_dens, zmin_dens);
            carrier_density.Spin_Down.Save_1D_Data("dens_1D_down.dat", dz_dens, zmin_dens);

            final_dens_solv.Output(carrier_density + dopent_density, "charge_density.dat", false);
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dz_pot, nz_pot, zmin_pot) - chem_pot, "potential.dat");

            double tot_dens_Quantum = (from val in carrier_density.Spin_Summed_Data.vec
                               where val < 0.0
                               select -1.0e14 * val * dz_dens / Physics_Base.q_e).ToArray().Sum();
            Console.WriteLine("Carrier density at heterostructure interface: \t" + tot_dens_Quantum.ToString("e3") + " cm^-2");
        }
    }
}