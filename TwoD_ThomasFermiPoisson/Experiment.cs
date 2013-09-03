using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    enum Density_Method
    {
        by_k,
        by_E,
        by_ThomasFermi
    }

    class Experiment : Experiment_Base
    {
        SpinResolved_DoubleMatrix charge_density;

        double dy, dz;
        int ny, nz; 

        Density_Method calculation_method = Density_Method.by_ThomasFermi;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // simulation domain inputs
            if (input_dict.ContainsKey("dy")) this.dy = (double)input_dict["dy"]; else throw new KeyNotFoundException("No dy in input dictionary!");
            if (input_dict.ContainsKey("ny")) this.ny = (int)(double)input_dict["ny"]; else throw new KeyNotFoundException("No ny in input dictionary!");
            if (input_dict.ContainsKey("dz")) this.dz = (double)input_dict["dz"]; else throw new KeyNotFoundException("No dz in input dictionary!");
            if (input_dict.ContainsKey("nz")) this.nz = (int)(double)input_dict["nz"]; else throw new KeyNotFoundException("No nz in input dictionary!");

            // physics parameters are done by the base method
            base.Initialise(input_dict, dz, nz);

            // try to get the potential and density from the dictionary... they probably won't be there and if not... make them
            if (input_dict.ContainsKey("SpinResolved_Density")) this.charge_density = (SpinResolved_DoubleMatrix)input_dict["SpinResolved_Density"]; else this.charge_density = new SpinResolved_DoubleMatrix(ny, nz);
            if (input_dict.ContainsKey("Potential")) this.band_offset = new Band_Data((DoubleMatrix)input_dict["Potential"]); else this.band_offset = new Band_Data(new DoubleMatrix(ny, nz));
        }

        public override void Run()
        {
            // create classes and initialise
            TwoD_ThomasFermiSolver dens_solv = new TwoD_ThomasFermiSolver(band_structure, acceptor_conc, donor_conc, acceptor_energy, donor_energy, temperature, dy, dz, ny, nz);
            double bottom_bc = dens_solv.Get_Chemical_Potential(nz - 1);

            TwoD_PoissonSolver pois_solv = new TwoD_PoissonSolver(dy, dz, ny, nz, bottom_bc, layer_depths, using_flexPDE, flexdPDE_input, freeze_dopents, tol);
            // initialise the band energy as the solution with zero density
            band_offset = new Band_Data(pois_solv.Get_Band_Energy(new DoubleMatrix(ny, nz, 0.0)));

            // run self consistent loop
            int count = 0;
            while (!pois_solv.Converged)
            {
                Console.WriteLine("Iteration:\t" + count.ToString());

                // calculate the total charge density for this band offset
                charge_density = dens_solv.Get_TwoD_ChargeDensity(band_offset.mat);

                // solve the band energy for the given charge density and mix in with the old band energy
                Band_Data new_band_energy = new Band_Data(pois_solv.Get_Band_Energy(charge_density.Spin_Summed_Matrix));
                pois_solv.Blend(ref band_offset, new_band_energy, alpha);

                // change the potential mixing parameter
                //if ((count + 1) % potential_mixing_rate == 0)
                //    alpha = pois_solv.Renew_Mixing_Parameter(potential, new_potential, alpha_prime, alpha);

                // transfer new potential array to potential array
                //potential = new_potential;
                count++;
            }

            dens_solv.Output(charge_density, "charge_density.dat");
            pois_solv.Output(Input_Band_Structure.Expand_BandStructure(band_structure, ny) / 2.0 - band_offset, "potential.dat");
        }
    }
}
