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
        SpinResolved_DoubleVector dopent_concentration;
        DoubleVector acceptor_concentration, donor_concentration;
        DoubleVector acceptor_energy_above_Ev, donor_energy_below_Ec;

        double no_kB_T_above_Ef = 10.0;
        double dE;

        public OneD_ThomasFermiSolver(DoubleVector band_gap, DoubleVector acceptor_concentration, DoubleVector donor_concentration, DoubleVector acceptor_energy, DoubleVector donor_energy,
                                        double fermi_Energy, double temperature, double dE, double dz, int nz) 
            : base(fermi_Energy, temperature, 1.0, 1.0, dz, 1, 1, nz)
        {
            this.dE = dE;

            // set band profile with spin-degeneracy
            this.band_gap = band_gap;

            // set dopent concentrations and copy to a vector of double the length for spin-degeneracy
            this.acceptor_concentration = acceptor_concentration;
            this.donor_concentration = donor_concentration;

            this.dopent_concentration = (SpinResolved_DoubleVector)(donor_concentration - acceptor_concentration);

            // and set relative dopent energies to the conduction band
            this.donor_energy_below_Ec = band_gap / 2.0 - donor_energy;
            this.acceptor_energy_above_Ev = band_gap / 2.0 - acceptor_energy;
        }

        public SpinResolved_DoubleVector Get_OneD_Density(DoubleVector conduction_band_energy)
        {
            // spin-resolved density
            SpinResolved_DoubleVector density = new SpinResolved_DoubleVector(nz);

            // Find conduction band density
            SpinResolved_DoubleVector conduction_density = Calculate_Conduction_Band_Density(conduction_band_energy);

            // Find valence band density
            DoubleVector valence_band = conduction_band_energy - band_gap;
            SpinResolved_DoubleVector valence_density = Calculate_Valence_Band_Density(valence_band);

            // Calculate donor occupation probability
            SpinResolved_DoubleVector acceptor_density = Calculate_Dopent_Density(- valence_band, acceptor_concentration, acceptor_energy_above_Ev);
            SpinResolved_DoubleVector donor_density = Calculate_Dopent_Density(conduction_band_energy, donor_concentration, donor_energy_below_Ec);

            // return total density (for 1D; hence the divide-by-lattice-spacing)
            return -1.0 * q_e * (dopent_concentration + valence_density + acceptor_density - conduction_density - donor_density);
        }

        /// <summary>
        /// calculates the conduction band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        SpinResolved_DoubleVector Calculate_Conduction_Band_Density(DoubleVector conduction_band_energy)
        {
            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(nz);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((no_kB_T_above_Ef * kB * temperature - conduction_band_energy.Min()) / dE);
            else
                no_energy_steps = 100;

            for (int j = 0; j < nz; j++)
                for (int i = 0; i < no_energy_steps; i++)
                {
                    double energy_above_eF = conduction_band_energy[j] + i * dE;

                    result[j, Spin.Up] += Get_3D_DensityofStates(conduction_band_energy[j], energy_above_eF) * Get_Fermi_Function(energy_above_eF) * dE;
                    result[j, Spin.Down] += Get_3D_DensityofStates(conduction_band_energy[j], energy_above_eF) * Get_Fermi_Function(energy_above_eF) * dE;
                }

            return result;
        }

        /// <summary>
        /// calculates the valence band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        SpinResolved_DoubleVector Calculate_Valence_Band_Density(DoubleVector valence_band)
        {
            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(nz);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((valence_band.Max() - no_kB_T_above_Ef * kB * temperature) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < no_energy_steps; i++)
                for (int j = 0; j < nz; j++)
                {
                    double energy_below_eF = valence_band[j] - i * dE;

                    // energies are negative (for density of states) as we are considering holes
                    // also, we subtract from "result" as holes are negatively charged
                    result[j, Spin.Up] += Get_3D_DensityofStates(-1.0 * valence_band[j], -1.0 * energy_below_eF) * (1.0 - Get_Fermi_Function(energy_below_eF)) * dE;
                    result[j, Spin.Down] += Get_3D_DensityofStates(-1.0 * valence_band[j], -1.0 * energy_below_eF) * (1.0 - Get_Fermi_Function(energy_below_eF)) * dE;
                }

            return result;
        }

        /// <summary>
        /// calculates the dopent density profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        SpinResolved_DoubleVector Calculate_Dopent_Density(DoubleVector conduction_band_energy, DoubleVector dopent_conc, DoubleVector dopent_energy_below_Ec)
        {
            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(nz);

            // calculates the occupation of the dopents 
            // (note that the dopent concentrations are only a function of depth)
            for (int i = 0; i < nz; i++)
            {
                result[i, Spin.Up] = dopent_conc[i] * Get_Dopent_Fermi_Function(conduction_band_energy[i] - dopent_energy_below_Ec[i]);
                result[i, Spin.Down] = dopent_conc[i] * Get_Dopent_Fermi_Function(conduction_band_energy[i] - dopent_energy_below_Ec[i]);
            }

            return result;
        }

        /// <summary>
        /// gets 3D spin-resolved density of states for given potential and energy
        /// </summary>
        double Get_3D_DensityofStates(double potential, double energy)
        {
            // calculate dk^2 / dE and k^2
            double dk2_dE = (2 * mass) / hbar * hbar;
            double k2 = (energy - potential) * dk2_dE;
            
            // prefactor for number of states per unit volume... ie (4 pi / 3) / (2 pi)^3
            double geometric_prefactor = 1.0 / (6 * Math.PI * Math.PI);

            double density_of_states = geometric_prefactor * dk2_dE * 1.5 * Math.Pow(k2, 0.5);
            return density_of_states;
        }

    }
}
