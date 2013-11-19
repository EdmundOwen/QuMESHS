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
        public void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

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

        public override void Run()
        {
            // get temperatures to run the experiment at
            double[] run_temps = Freeze_Out_Temperatures();
            
            for (int i = 0; i < run_temps.Length; i++)
            {
                double current_temperature = run_temps[i];

                // create charge density solver and calculate boundary conditions
                OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(current_temperature, dz, nz, zmin);
                double top_bc = 0;// dens_solv.Get_Chemical_Potential(0);
                double bottom_bc = dens_solv.Get_Chemical_Potential(layers, zmin);

                // initialise potential solver
                OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, layers, using_flexPDE, flexPDE_input, flexPDE_location, tol);
                pois_solv.Set_Boundary_Conditions(layers, top_bc, bottom_bc);

                // initialise the band energy as the solution with zero density
                band_offset = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data);

                int count = 0;
                while (!pois_solv.Converged)
                {
                    Console.WriteLine("Iteration: " + count.ToString() + "\ttemperature: " + current_temperature.ToString() + "\tConvergence factor: " + pois_solv.Convergence_Factor.ToString());

                    // calculate the total charge density for this band offset
                    dens_solv.Get_OneD_ChargeDensity(layers, ref charge_density, band_offset);

                    // solve the band energy for the givencharge  density and mix in with the old band energy
                    Band_Data new_band_energy = pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Data);
                    pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                    count++;
                }
            }

            // initialise output solvers
            OneD_ThomasFermiSolver final_dens_solv = new OneD_ThomasFermiSolver(temperature, dz, nz, zmin);
            OneD_PoissonSolver final_pois_solv = new OneD_PoissonSolver(dz, nz, layers, using_flexPDE, flexPDE_input, flexPDE_location, tol);

            // save final density out
            final_pois_solv.Save_Density(charge_density.Spin_Summed_Data, "dens_1D.dat");

            final_dens_solv.Output(charge_density, "charge_density.dat");
            final_pois_solv.Output(Input_Band_Structure.Get_BandStructure_Grid(layers, dz, nz, zmin) - band_offset, "potential.dat");
        }
    }
}