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
    class Experiment
    {
        bool converged = false;
        double alpha, alpha_prime, tol;
        int potential_mixing_rate = 10;
        
        bool using_flexPDE = false;
        string flexdPDE_input;

        Potential_Data potential;
        SpinResolved_DoubleVector density;
        DoubleVector band_structure;
        DoubleVector acceptor_energy, donor_energy;
        DoubleVector acceptor_conc, donor_conc;

        double temperature = 300.0;

        double dz;
        int nz;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

            // solver inputs
            if (input_dict.ContainsKey("tolerance")) this.tol = (double)input_dict["tolerance"]; else throw new KeyNotFoundException("No solution tolerance in input dictionary!");
            if (input_dict.ContainsKey("alpha")) { this.alpha = (double)input_dict["alpha"]; alpha_prime = alpha; } else throw new KeyNotFoundException("No potential mixing parameter, alpha, in input dictionary!");
            if (input_dict.ContainsKey("FlexPDE_file")) { this.using_flexPDE = true; this.flexdPDE_input = (string)input_dict["FlexPDE_file"]; }

            // physical inputs
            if (input_dict.ContainsKey("T")) this.temperature = (double)input_dict["T"]; else throw new KeyNotFoundException("No temperature in input dictionary!");

            // get the band structure
            if (input_dict.ContainsKey("BandStructure_File"))
            {
                band_structure = Input_Band_Structure.GetBandStructure((string)input_dict["BandStructure_File"], nz, dz);

                Input_Band_Structure.GetDopentData((string)input_dict["BandStructure_File"], nz, dz, Dopent.acceptor, out acceptor_conc, out acceptor_energy);
                Input_Band_Structure.GetDopentData((string)input_dict["BandStructure_File"], nz, dz, Dopent.donor, out donor_conc, out donor_energy);
            }
            else throw new KeyNotFoundException("No band structure file found in input dictionary!");

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.density = (SpinResolved_DoubleVector)input_dict["SpinResolved_Density"]; else this.density = new SpinResolved_DoubleVector(nz);
            if (input_dict.ContainsKey("Potential")) this.potential = new Potential_Data((DoubleMatrix)input_dict["Potential"]); else potential = new Potential_Data(new DoubleVector(nz));
        }

        public void Initialise(DoubleVector Potential, DoubleVector Density, double dz, double alpha, double tol, int nz)
        {
            throw new NotImplementedException();
            //this.potential = Potential; this.density = Density;
            Initialise(dz, alpha, tol, nz);
        }

        public void Initialise(double dz, double alpha, double tol, int nz)
        {
            this.alpha = alpha; this.alpha_prime = alpha;

            this.tol = tol;

            this.nz = nz; this.dz = dz;
        }

        public void Run()
        {
            // create density solver and calculate boundary conditions
            OneD_ThomasFermiSolver dens_solv = new OneD_ThomasFermiSolver(band_structure, acceptor_conc, donor_conc, acceptor_energy, donor_energy, 0.0, temperature, 10.0, dz, nz);
            double top_bc = dens_solv.Get_Chemical_Potential(0);
            double bottom_bc = dens_solv.Get_Chemical_Potential(nz - 1);

            // initialise potential solver
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, top_bc, bottom_bc, using_flexPDE, flexdPDE_input, tol);

            int count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString());

                // calculate the total density for this potential
                density = dens_solv.Get_OneD_Density(band_structure / 2.0 + potential.vec);

                // solve the potential for the given density and mix in with the old potential
                Potential_Data new_potential = pois_solv.Blend(potential, new Potential_Data(pois_solv.Get_Potential(density.Spin_Summed_Vector)), alpha);

                // change the potential mixing parameter
                if ((count + 1) % potential_mixing_rate == 0)
                    alpha = pois_solv.Renew_Mixing_Parameter(potential, new_potential, alpha_prime, alpha);
                    //alpha = pois_solv.Renew_Mixing_Factor(potential, new_potential);

                // transfer new potential array to potential array
                potential = new_potential;
                count++;
            }

            Solver_Bases.ZeroD_Charge tmp = new ZeroD_Charge(band_structure[0], acceptor_conc[0], acceptor_energy[0], donor_conc[0], donor_energy[0], temperature);
            
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

            dens_solv.Output(density, "density.dat");
            pois_solv.Output(new Potential_Data(potential.vec + band_structure / 2.0), "potential.dat");
        }
    }
}