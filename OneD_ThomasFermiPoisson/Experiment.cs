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
        double alpha, alpha_prime, tol;
        
        bool using_flexPDE = false;
        string flexdPDE_input;

        Band_Data band_offset;
        SpinResolved_DoubleVector density;

        double[] layer_depths;
        DoubleVector band_structure;
        DoubleVector acceptor_energy, donor_energy;
        DoubleVector acceptor_conc, donor_conc;

        double temperature = 300.0;
        
        // temperature at which the dopents are assumed to have frozen out
        double freeze_out_T = 70.0;
        bool freeze_dopents;

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
            // and check whether the dopents are frozen out
            if (input_dict.ContainsKey("dopents_frozen")) this.freeze_dopents = (bool)input_dict["dopents_frozen"]; 
            else { if (temperature < freeze_out_T) freeze_dopents = true; else freeze_dopents = false;}

            // get the band structure
            if (input_dict.ContainsKey("BandStructure_File"))
            {
                Input_Band_Structure band_structure_generator = new Input_Band_Structure((string)input_dict["BandStructure_File"]);

                band_structure = band_structure_generator.GetBandStructure(nz, dz);
                layer_depths = band_structure_generator.Layer_Depths;

                band_structure_generator.GetDopentData(nz, dz, Dopent.acceptor, out acceptor_conc, out acceptor_energy);
                band_structure_generator.GetDopentData(nz, dz, Dopent.donor, out donor_conc, out donor_energy);
            }
            else throw new KeyNotFoundException("No band structure file found in input dictionary!");

            // try to get the band offset and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.density = (SpinResolved_DoubleVector)input_dict["SpinResolved_Density"]; else this.density = new SpinResolved_DoubleVector(nz);
            if (input_dict.ContainsKey("Band_Offset")) this.band_offset = new Band_Data((DoubleMatrix)input_dict["Potential"]); else band_offset = new Band_Data(new DoubleVector(nz));
        }

        public void Initialise(DoubleVector Band_Offset, DoubleVector Density, double dz, double alpha, double tol, int nz)
        {
            throw new NotImplementedException();
            //this.Band_Offset = Band_Offset; this.density = Density;
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
            double top_bc = 0;// dens_solv.Get_Chemical_Potential(0);
            double bottom_bc = 700;// dens_solv.Get_Chemical_Potential(nz - 1);

            // initialise potential solver
            OneD_PoissonSolver pois_solv = new OneD_PoissonSolver(dz, nz, top_bc, bottom_bc, layer_depths, using_flexPDE, flexdPDE_input, freeze_dopents, tol);

            int count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString());

                // calculate the total density for this band offset
                density = dens_solv.Get_OneD_Density(band_offset.vec);

                // solve the band energy for the given density and mix in with the old band energy
                Band_Data new_band_energy = new Band_Data(pois_solv.Get_Band_Energy(density.Spin_Summed_Vector));
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

            dens_solv.Output(density, "density.dat");
            pois_solv.Output(new Band_Data(band_structure / 2.0 - band_offset.vec), "potential.dat");
        }
    }
}