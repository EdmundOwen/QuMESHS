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
        int nz_dft;
        double zmin_dft, dz_dft;

        public new void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

            // dft inputs
            if (input_dict.ContainsKey("dft"))
                if (bool.Parse((string)input_dict["dft"]))
                {
                    if (input_dict.ContainsKey("dz_dft")) this.dz_dft = (double)input_dict["dz_dft"]; else throw new KeyNotFoundException("No dz for dft calculation in input dictionary!");
                    if (input_dict.ContainsKey("nz_dft")) this.nz_dft = (int)(double)input_dict["nz_dft"]; else throw new KeyNotFoundException("No nz for dft calculation in input dictionary!");
                    if (input_dict.ContainsKey("zmin_dft")) this.zmin_dft = (double)input_dict["zmin_dft"]; else throw new KeyNotFoundException("No minimum z value for dft calculation in input dictionary!");
                }

            // physics parameters are done by the base method
            base.Initialise(input_dict);

            // try to get the band offset and the charge density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.charge_density = (SpinResolved_Data)input_dict["SpinResolved_Density"];
            else this.charge_density = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));

            if (input_dict.ContainsKey("Band_Offset")) this.band_offset = new Band_Data((DoubleVector)input_dict["Potential"]); else band_offset = new Band_Data(new DoubleVector(nz));
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
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, layers, using_flexPDE, flexPDE_input, flexPDE_location, tol);
            
            for (int i = 0; i < run_temps.Length; i++)
            {
                double current_temperature = run_temps[i];

                // create charge density solver and calculate boundary conditions
                OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(current_temperature, dz, nz, zmin);
                double top_bc = 0.0;// dens_solv.Get_Chemical_Potential(0);
                double bottom_bc = dens_solv.Get_Chemical_Potential(zmin, layers);

                // set the boundary conditions
                pois_solv.Set_Boundary_Conditions(layers, top_bc, bottom_bc);

                // initialise the band energy as the solution with zero density
                band_offset = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data);

                count = 0;
                while (!pois_solv.Converged)
                {
                    Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + current_temperature.ToString() + "\tConvergence factor: " + pois_solv.Convergence_Factor.ToString());

                    // calculate the total charge density for this band offset
                    dens_solv.Get_ChargeDensity(layers, ref charge_density, band_offset);

                    // solve the band energy for the givencharge  density and mix in with the old band energy
                    Band_Data new_band_energy = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data);
                    pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                    count++;
                }

                pois_solv.Reset();
            }

            // and then run the DFT solver at the base temperature over a limited range
            OneD_DFTSolver dft_solv = new OneD_DFTSolver(temperature, dz_dft, nz_dft, zmin_dft);
            SpinResolved_Data dft_dens = new SpinResolved_Data(new Band_Data(new DoubleVector(nz_dft)), new Band_Data(new DoubleVector(nz_dft)));
            // but with a smaller mixing parameter
            Console.WriteLine("Renewing mixing parameter for DFT calculation");
            alpha /= 10.0;


            count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + run_temps[run_temps.Length - 1].ToString() + "\tConvergence factor: " + pois_solv.Convergence_Factor.ToString());

                // calculate the total charge density for this band offset
                Band_Data dft_band_offset = new Band_Data(new DoubleVector(nz_dft));
                Interpolate_DFT_Grid(ref dft_dens, ref dft_band_offset, charge_density, band_offset);
                Get_Potential(ref dft_band_offset, layers);
                dft_solv.Get_ChargeDensity(layers, ref dft_dens, dft_band_offset);

                // find the top and bottom of the dft solution domain for charge_density
                int offset_min = (int)Math.Round((zmin_dft - zmin) / dz);
                int offset_max = (int)Math.Round((zmin_dft + nz_dft * dz_dft - zmin) / dz);

                // and insert the dft density into charge_density
                for (int i = offset_min; i < offset_max; i++)
                {
                    // calculate the index in the dft array for the spatial positions of charge_density
                    int dft_index;
                    if (i == offset_max - 1)
                        dft_index = (int)Math.Floor((i - offset_min) * dz / dz_dft);
                    else
                        dft_index = (int)Math.Round((i - offset_min) * dz / dz_dft);

                    charge_density.Spin_Up[i] = dft_dens.Spin_Up[dft_index];
                    charge_density.Spin_Down[i] = dft_dens.Spin_Down[dft_index];
                }

                // solve the band energy for the given charge_density and mix in with the old band energy
                Band_Data new_band_energy = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data);
                pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                count++;
            }

            pois_solv.Reset();

            // initialise output solvers
            OneD_ThomasFermiSolver final_dens_solv = new OneD_ThomasFermiSolver(temperature, dz, nz, zmin);
            OneD_PoissonSolver final_pois_solv = new OneD_PoissonSolver(dz, nz, layers, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            final_pois_solv.Save_Density(charge_density.Spin_Summed_Data, "dens_1D.dat");

            final_dens_solv.Output(charge_density, "charge_density.dat");
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dz, nz, zmin) - band_offset, "potential.dat");
        }

        void Get_Potential(ref Band_Data dft_band_offset, ILayer[] layers)
        {
            for (int i = 0; i < nz_dft; i++)
            {
                double pos = zmin_dft + i * dz_dft;
                double band_gap = Geom_Tool.GetLayer(layers, pos).Band_Gap;
                dft_band_offset[i] = 0.5 * band_gap - dft_band_offset[i];
            }
        }

        void Interpolate_DFT_Grid(ref SpinResolved_Data dft_dens, ref Band_Data dft_band_offset, SpinResolved_Data charge_density, Band_Data band_offset)
        {
            for (int i = 0; i < nz_dft; i++)
            {
                double dft_pos = zmin_dft + i * dz_dft;
                int init_index = (int)Math.Floor((i * dz_dft - zmin + zmin_dft) / dz); // index for the initial domain (ie. for charge_density and band_offset)

                dft_dens.Spin_Up[i] = charge_density.Spin_Up[init_index] + (charge_density.Spin_Up[init_index + 1] - charge_density.Spin_Up[init_index]) * (dft_pos - Math.Floor(dft_pos / dz) * dz) / dz;
                dft_dens.Spin_Down[i] = charge_density.Spin_Down[init_index] + (charge_density.Spin_Down[init_index + 1] - charge_density.Spin_Down[init_index]) * (dft_pos - Math.Floor(dft_pos / dz) * dz) / dz;
                dft_band_offset[i] = band_offset[init_index] + (band_offset[init_index + 1] - band_offset[init_index]) * (dft_pos - Math.Floor(dft_pos / dz) * dz) / dz;
            }
        }
    }
}