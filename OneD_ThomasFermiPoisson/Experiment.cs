using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Analysis;
using Solver_Bases;

namespace OneD_ThomasFermiPoisson
{
    class Experiment : Experiment_Base
    {
        SpinResolved_DoubleVector charge_density;

        double dz;
        int nz;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

            // physics parameters are done by the base method
            base.Initialise(input_dict, dz, nz);

            // try to get the band offset and the charge density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.charge_density = (SpinResolved_DoubleVector)input_dict["SpinResolved_Density"]; else this.charge_density = new SpinResolved_DoubleVector(nz);
            if (input_dict.ContainsKey("Band_Offset")) this.band_offset = new Band_Data((DoubleMatrix)input_dict["Potential"]); else band_offset = new Band_Data(new DoubleVector(nz));
        }

        public override void Run()
        {
            // create charge density solver and calculate boundary conditions
            OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(band_structure, acceptor_conc, donor_conc, acceptor_energy, donor_energy, temperature, dz, nz);
            double top_bc = 0;// dens_solv.Get_Chemical_Potential(0);
            double bottom_bc = dens_solv.Get_Chemical_Potential(nz - 1);

            // initialise potential solver
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, top_bc, bottom_bc, layer_depths, using_flexPDE, flexdPDE_input, freeze_dopents, tol);
            // initialise the band energy as the solution with zero density
            band_offset = new Band_Data(pois_solv.Get_Band_Energy(new DoubleVector(nz, 0.0)));

            int count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString());

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
            */

            dens_solv.Output(charge_density, "charge_density.dat");
            pois_solv.Output(new Band_Data(band_structure / 2.0 - band_offset.vec), "potential.dat");
        }
    }
}