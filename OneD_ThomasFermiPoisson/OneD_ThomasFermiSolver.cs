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
        DoubleVector acceptor_energy_below_Ec, donor_energy_below_Ec;

        double no_kB_T_above_Ef = 100.0;
        double dE;

        public OneD_ThomasFermiSolver(DoubleVector band_gap, DoubleVector acceptor_concentration, DoubleVector donor_concentration, DoubleVector acceptor_energy, DoubleVector donor_energy,
                                        double fermi_Energy, double temperature, double dE, int nz) 
            : base(fermi_Energy, temperature, 1, 1, nz)
        {
            this.dE = dE;

            // set band profile with spin-degeneracy
            this.band_gap = band_gap;

            // set dopent concentrations and copy to a vector of double the length for spin-degeneracy
            this.acceptor_concentration = new DoubleVector(2 * nz);
            this.donor_concentration = new DoubleVector(2 * nz);
            for (int i = 0; i < nz; i++)
            {
                this.acceptor_concentration[i] = acceptor_concentration[i]; this.acceptor_concentration[i + nz] = acceptor_concentration[i];
                this.donor_concentration[i] = donor_concentration[i]; this.donor_concentration[i + nz] = donor_concentration[i];
            }

            // and set relative dopent energies to the conduction band
            this.donor_energy_below_Ec = band_gap / 2.0 - donor_energy;
            this.acceptor_energy_below_Ec = band_gap / 2.0 - acceptor_energy;
        }

        public DoubleVector Get_OneD_Density(DoubleVector conduction_band_energy)
        {
            // spin-resolved density
            DoubleVector density = new DoubleVector(2 * nz);

            // Find conduction band density
            DoubleVector conduction_density = Calculate_Conduction_Band_Density(conduction_band_energy);

            // Find valence band density
            DoubleVector valence_density = Calculate_Valence_Band_Density(conduction_band_energy);

            // Calculate donor occupation probability
            DoubleVector acceptor_density = Calculate_Dopent_Density(conduction_band_energy, acceptor_concentration, acceptor_energy_below_Ec);
            DoubleVector donor_density = Calculate_Dopent_Density(conduction_band_energy, donor_concentration, donor_energy_below_Ec);

            // return total density
            return -1.0 *  q_e * (donor_concentration - acceptor_concentration + valence_density + acceptor_density - conduction_density - donor_density);
        }

        /// <summary>
        /// calculates the conduction band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        DoubleVector Calculate_Conduction_Band_Density(DoubleVector conduction_band_energy)
        {
            DoubleVector result = new DoubleVector(2 * nz);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((no_kB_T_above_Ef * kB * temperature - conduction_band_energy.Min()) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < no_energy_steps; i++)
                for (int j = 0; j < 2 * nz; j++)
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
            DoubleVector result = new DoubleVector(2 * nz);
            DoubleVector valence_band = conduction_band_energy - band_gap;

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((valence_band.Max() - no_kB_T_above_Ef * kB * temperature) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < no_energy_steps; i++)
                for (int j = 0; j < 2 * nz; j++)
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
        DoubleVector Calculate_Dopent_Density(DoubleVector conduction_band_energy, DoubleVector dopent_concentration, DoubleVector dopent_energy_below_Ec)
        {
            DoubleVector result = new DoubleVector(2 * nz);

            // calculate today density of dopents (dopent charge and carrier charge)
            for (int i = 0; i < nz; i++)
            {
                // spin-up
                result[i] = dopent_concentration[i] * Get_Dopent_Occupation(conduction_band_energy[i] - dopent_energy_below_Ec[i]);
                // and spin-down is degenerate here
                result[i + nz] = dopent_concentration[i] * Get_Dopent_Occupation(conduction_band_energy[i] - dopent_energy_below_Ec[i]);
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
