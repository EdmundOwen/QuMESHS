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

            // see whether there's a top gate voltage for determining top_bc
            Get_From_Dictionary<double>(input_dict, "top_V", ref top_V, true);

            // try to get the chemical potential and the charge density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.charge_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
            else this.charge_density = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_pot)), new Band_Data(new DoubleVector(nz_pot)));

            if (input_dict.ContainsKey("Chemical_Potential")) this.chem_pot = new Band_Data((DoubleVector)input_dict["Chemical_Potential"]); else chem_pot = new Band_Data(new DoubleVector(nz_pot));

            if (input_dict.ContainsKey("dft")) this.TF_only = !bool.Parse((string)input_dict["dft"]);
            if (input_dict.ContainsKey("TF_only")) this.TF_only = bool.Parse((string)input_dict["TF_only"]);
        }

        /*
        public override void Run()
        {
            // create charge density solver and calculate boundary conditions
            OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(band_structure, acceptor_conc, donor_conc, acceptor_energy, donor_energy, temperature, dz, nz);
            double top_bc = 0;// dens_solv.Get_Chemical_Potential(0);
            double bottom_bc = dens_solv.Get_Chemical_Potential(nz - 1);

            // initialise potential solver
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, top_bc, bottom_bc, layer_depths, using_flexPDE, flexPDE_input, flexPDE_location, freeze_dopents, tol);
            // initialise the band energy as the solution with zero density
            band_offset = new Band_Data(pois_solv.Get_Band_Energy(new DoubleVector(nz, 0.0)));

            int count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString() + "\tConvergence factor:\t" + pois_solv.Convergence_Factor.ToString());

                // calculate the total charge density for this band offset
                charge_density = dens_solv.Get_OneD_ChargeDensity(band_offset.vec);

                // solve the band energy for the givencharge  density and mix in with the old band energy
                Band_Data new_band_energy = new Band_Data(pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Vector));
                pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                // change the band energy mixing parameter
                //if ((count + 1) % potential_mixing_rate == 0)
                //    alpha = pois_solv.Renew_Mixing_Parameter(band_offset, new_band_energy, alpha_prime, alpha);
                    //alpha = pois_solv.Renew_Mixing_Factor(band_offset, new_band_energy);

                // transfer new band energy array to band energy array
                //band_offset = new_band_energy;
                count++;
            }

            /*
            OneD_Error_Functional error_func = new OneD_Error_Functional(dens_solv, pois_solv, nz);
            List<Constraint> constraints = new List<Constraint>();
            DoubleFunctionalDelegate c1 =            new DoubleFunctionalDelegate(nz, new Func<DoubleVector,double>(delegate(DoubleVector v) { return 0.0; }));
            NonlinearConstraint constraint1 = new NonlinearConstraint(            c1, ConstraintType.GreaterThanOrEqualTo);
            constraints.Add(constraint1);
            NonlinearProgrammingProblem problem = new NonlinearProgrammingProblem(error_func, constraints);

            ActiveSetLineSearchSQP solver = new ActiveSetLineSearchSQP(20000, tol);
            converged = solver.Solve(problem, new DoubleVector(nz));
            Console.WriteLine("Termination status = " + solver.SolverTerminationStatus);
            density.Spin_Up = solver.OptimalX;
            potential.vec = pois_solv.Get_Potential(density.Spin_Summed_Vector);
            

            // save final density out
            pois_solv.Save_Density(new Band_Data(charge_density.Spin_Summed_Vector), "dens_1D.dat");

            dens_solv.Output(charge_density, "charge_density.dat");
            pois_solv.Output(new Band_Data(band_structure / 2.0 - band_offset.vec), "potential.dat");
        }*/

        int count = 0;
        public override void Run()
        {
            // get temperatures to run the experiment at
            double[] run_temps = Freeze_Out_Temperatures();

            // initialise potential solver
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            
            for (int i = 0; i < run_temps.Length; i++)
            {
                double current_temperature = run_temps[i];

                // create charge density solver and calculate boundary conditions
                OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(current_temperature, dz_pot, nz_pot, zmin_pot);
                double top_bc = 0.0;// dens_solv.Get_Chemical_Potential(0);
                this.bottom_bc = dens_solv.Get_Chemical_Potential(zmin_pot, layers);

                // set the boundary conditions
                pois_solv.Set_Boundary_Conditions(layers, top_bc, bottom_bc, Geom_Tool.Get_Zmin(layers) + dz_pot * nz_pot, Geom_Tool.Get_Zmin(layers));

                // initialise the chemical potential as the solution with zero density
                chem_pot = pois_solv.Get_Chemical_Potential(charge_density.Spin_Summed_Data);

                count = 0;
                while (!pois_solv.Converged)
                {
                    if (count % 100 == 0)
                        Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + current_temperature.ToString() + "\tConvergence factor: " + pois_solv.Convergence_Factor.ToString());

                    // calculate the total charge density for this chemical_potential
                    dens_solv.Get_ChargeDensity(layers, ref charge_density, chem_pot);

                    // solve the potential for the givencharge  density and mix in with the old chemical potential
                    Band_Data new_chem_pot = pois_solv.Get_Chemical_Potential(charge_density.Spin_Summed_Data);
                    pois_solv.Blend(ref chem_pot, new_chem_pot, alpha);

                    count++;
                }

                pois_solv.Reset();
            }

            if (TF_only)
                return;

            // and then run the DFT solver at the base temperature over a limited range
            OneD_DFTSolver dft_solv = new OneD_DFTSolver(temperature, dz_dens, nz_dens, zmin_dens);
            SpinResolved_Data dft_dens = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_dens)), new Band_Data(new DoubleVector(nz_dens)));
            // but with a smaller mixing parameter
            alpha /= 3.0;
            Console.WriteLine("Starting DFT calculation");

            // reset boundary conditions
            double tmp_bc = top_V * (Physics_Base.q_e * Physics_Base.energy_V_to_meVpzC);
            pois_solv.Set_Boundary_Conditions(layers, tmp_bc, bottom_bc, Geom_Tool.Get_Zmin(layers) + dz_pot * nz_pot, Geom_Tool.Get_Zmin(layers));

            count = 0;
            while (!pois_solv.Converged)
            {
                if (count % 100 == 0)
                    Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + run_temps[run_temps.Length - 1].ToString() + "\tConvergence factor: " + pois_solv.Convergence_Factor.ToString());

                // calculate the total charge density for this chemical potential
                Band_Data dft_chem_pot = new Band_Data(new DoubleVector(nz_dens));
                Interpolate_DFT_Grid(ref dft_dens, ref dft_chem_pot, charge_density, chem_pot);
                Get_Potential(ref dft_chem_pot, layers);
                dft_solv.Get_ChargeDensity(layers, ref dft_dens, dft_chem_pot);

                // find the top and bottom of the dft solution domain for charge_density
                int offset_min = (int)Math.Round((zmin_dens - zmin_pot) / dz_pot);
                int offset_max = (int)Math.Round((zmin_dens + nz_dens * dz_dens - zmin_pot) / dz_pot);

                // and insert the dft density into charge_density
                for (int i = offset_min; i < offset_max; i++)
                {
                    // calculate the index in the dft array for the spatial positions of charge_density
                    int dft_index;
                    if (i == offset_max - 1)
                        dft_index = (int)Math.Floor((i - offset_min) * dz_pot / dz_dens);
                    else
                        dft_index = (int)Math.Round((i - offset_min) * dz_pot / dz_dens);

                    charge_density.Spin_Up[i] = dft_dens.Spin_Up[dft_index];
                    charge_density.Spin_Down[i] = dft_dens.Spin_Down[dft_index];
                }

                // solve the chemical potential for the given charge_density and mix in with the old chemical potential
                Band_Data new_chem_pot = pois_solv.Get_Chemical_Potential(charge_density.Spin_Summed_Data);
                pois_solv.Blend(ref chem_pot, new_chem_pot, alpha);

                count++;
            }

            /*
            // Calculate the required surface charge needed to screen the dopents.
            // This will generate zero electric field above the structure
            double surface_charge = pois_solv.Get_Surface_Charge(band_offset, layers);
            // generate a new charge distribution
            SpinResolved_Data charge_w_surface = new SpinResolved_Data(new Band_Data(new DoubleVector(141)), new Band_Data(new DoubleVector(141)));
            for (int i = 1; i < 101; i++)
            {
                charge_w_surface.Spin_Up[i] = charge_density.Spin_Up[i];
                charge_w_surface.Spin_Down[i] = charge_density.Spin_Down[i];
            }
            charge_w_surface.Spin_Down[100] = surface_charge;

            pois_solv.Reset();

            // initialise output solvers
            OneD_ThomasFermiSolver final_dens_solv = new OneD_ThomasFermiSolver(temperature, dz, nz, zmin);
            OneD_PoissonSolver final_pois_solv = new OneD_PoissonSolver(dz, nz + 40, layers, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            Band_Data new_band_offset = final_pois_solv.Get_Band_Energy(charge_w_surface.Spin_Summed_Data);
           final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dz, nz+40, zmin) - new_band_offset, "pot_wsurface.dat");
             * */

            pois_solv.Reset();

            // initialise output solvers
            OneD_ThomasFermiSolver final_dens_solv = new OneD_ThomasFermiSolver(temperature, dz_dens, nz_dens, zmin_dens);
            OneD_PoissonSolver final_pois_solv = new OneD_PoissonSolver(this, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            //charge_density.Spin_Summed_Data.Save_1D_Data("dens_1D.dat", dz_dens, zmin_dens);

            //StreamWriter sw = new StreamWriter(dft_solv.No_Wavefunctions.ToString());
            //sw.WriteLine("Number of wave functions calculated = " + dft_solv.No_Wavefunctions.ToString());
            //sw.Close();
            //final_dens_solv.Output(charge_density / Physics_Base.q_e, "density.dat", false);

            final_dens_solv.Output(charge_density, "charge_density.dat", false);
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dz_pot, nz_pot, zmin_pot) - chem_pot, "potential.dat");
        }

        void Get_Potential(ref Band_Data dft_band_offset, ILayer[] layers)
        {
            for (int i = 0; i < nz_dens; i++)
            {
                double pos = zmin_dens + i * dz_dens;
                double band_gap = Geom_Tool.GetLayer(layers, pos).Band_Gap;
                dft_band_offset[i] = 0.5 * band_gap - dft_band_offset[i];
            }
        }

        void Interpolate_DFT_Grid(ref SpinResolved_Data dft_dens, ref Band_Data dft_band_offset, SpinResolved_Data charge_density, Band_Data band_offset)
        {
            for (int i = 0; i < nz_dens; i++)
            {
                double dft_pos = zmin_dens + i * dz_dens;
                int init_index = (int)Math.Floor((i * dz_dens - zmin_pot + zmin_dens) / dz_pot); // index for the initial domain (ie. for charge_density and band_offset)

                // no interpolation (it doesn't work...)
                dft_dens.Spin_Up[i] = charge_density.Spin_Up[init_index];// +(charge_density.Spin_Up[init_index + 1] - charge_density.Spin_Up[init_index]) * (dft_pos - Math.Floor(dft_pos / dz_pot) * dz_pot) / dz_pot;
                dft_dens.Spin_Down[i] = charge_density.Spin_Down[init_index];// +(charge_density.Spin_Down[init_index + 1] - charge_density.Spin_Down[init_index]) * (dft_pos - Math.Floor(dft_pos / dz_pot) * dz_pot) / dz_pot;
                dft_band_offset[i] = band_offset[init_index];// +(band_offset[init_index + 1] - band_offset[init_index]) * (dft_pos - Math.Floor(dft_pos / dz_pot) * dz_pot) / dz_pot;
            }
        }
    }
}