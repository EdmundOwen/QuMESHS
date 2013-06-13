using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace OneD_ThomasFermiPoisson
{
    class OneD_ThomasFermiSolver : Density_Solver
    {
        DoubleVector band_gap;
        DoubleVector acceptor_concentration, donor_concentration;
        double acceptor_energy_below_Ec, donoe_energy_below_Ec;

        double no_kB_T_above_Ef = 10.0;
        double dE;

        public OneD_ThomasFermiSolver(double fermi_Energy, double temperature, int nx, int ny, int nz) 
            : base(fermi_Energy, temperature, nx, ny, nz)
        {
        }

        public DoubleVector Get_OneD_Density(DoubleVector conduction_band_energy)
        {
            // spin-resolved density
            DoubleVector density = new DoubleVector(2 * nz);

            // Find conduction band density
            density = Calculate_Conduction_Band_Density(conduction_band_energy);

            // Find valence band density
            density += Calculate_Valence_Band_Density(conduction_band_energy);

            // Calculate donor occupation probability
            density -= Calculate_Dopent_Density(conduction_band_energy, acceptor_concentration, acceptor_energy_below_Ec);
            density += Calculate_Dopent_Density(conduction_band_energy, donor_concentration, donoe_energy_below_Ec);

            // return total density
            return density;
        }

        /// <summary>
        /// calculates the conduction band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        DoubleVector Calculate_Conduction_Band_Density(DoubleVector conduction_band_energy)
        {
            DoubleVector result = new DoubleVector(2 * nx);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((no_kB_T_above_Ef * kB * temperature - conduction_band_energy.Min()) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < no_energy_steps; i++)
                for (int j = 0; j < 2 * nx; j++)
                {
                    double energy_above_eF = conduction_band_energy[j % nz] + i * dE;

                    result[j] += Get_3D_DensityofStates(conduction_band_energy[j % nz], energy_above_eF) * Get_Fermi_Function(energy_above_eF);
                }

            return result;
        }

        /// <summary>
        /// calculates the valence band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        DoubleVector Calculate_Valence_Band_Density(DoubleVector conduction_band_energy)
        {
            DoubleVector result = new DoubleVector(2 * nx);
            DoubleVector valence_band = conduction_band_energy - band_gap;

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)(-1.0 * (valence_band.Max() - no_kB_T_above_Ef * kB * temperature) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < no_energy_steps; i++)
                for (int j = 0; j < 2 * nx; j++)
                {
                    double energy_below_eF = valence_band[j % nz] - i * dE;

                    // energies are negative (for density of states) as we are considering holes
                    // also, we subtract from "result" as holes are negatively charged
                    result[j] -= Get_3D_DensityofStates(-1.0 * valence_band[j % nz], -1.0 * energy_below_eF) * (1.0 - Get_Fermi_Function(energy_below_eF));
                }

            return result;
        }

        /// <summary>
        /// calculates the dopent density profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        DoubleVector Calculate_Dopent_Density(DoubleVector conduction_band_energy, DoubleVector dopent_concentration, double dopent_energy_below_Ec)
        {
            DoubleVector result = new DoubleVector(2 * nz);

            // calculate today density of dopents (dopent charge and carrier charge)
            for (int i = 0; i < nz; i++)
            {
                // spin-up
                result[i] = dopent_concentration[i] * (1.0 - Get_Dopent_Occupation(conduction_band_energy[i] - dopent_energy_below_Ec));
                // and spin-down is degenerate here
                result[i + nz] = dopent_concentration[i] * (1.0 - Get_Dopent_Occupation(conduction_band_energy[i] - dopent_energy_below_Ec));
            }

            return result;
        }

        /// <summary>
        /// gets 3D spin-resolved density of states for given potential and energy
        /// </summary>
        double Get_3D_DensityofStates(double potential, double energy)
        {
            // calculate dk^2 / dE and k^2
            double dk2_dE = hbar * hbar / (2 * mass);
            double k2 = (energy - potential) * dk2_dE;
            
            // prefactor for number of states per unit volume... ie (4 pi / 3) / (2 pi)^3
            double geometric_prefactor = 1.0 / (6 * Math.PI * Math.PI);

            double density_of_states = geometric_prefactor * dk2_dE * 1.5 * Math.Pow(k2, 0.5);
            return density_of_states;
        }

    }
}
