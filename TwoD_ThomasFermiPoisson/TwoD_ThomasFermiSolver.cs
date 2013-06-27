using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_ThomasFermiSolver : Density_Solver
    {
        DoubleMatrix band_gap;
        SpinResolved_DoubleMatrix dopent_concentration;
        DoubleVector acceptor_concentration, donor_concentration;
        DoubleVector acceptor_energy_below_Ec, donor_energy_below_Ec;

        double no_kB_T_above_Ef = 100.0;
        double dE;

        public TwoD_ThomasFermiSolver(DoubleVector band_gap, DoubleVector acceptor_concentration, DoubleVector donor_concentration, DoubleVector acceptor_energy, DoubleVector donor_energy,
                                        double fermi_Energy, double temperature, double dE, double dy, double dz, int ny, int nz) 
            : base(fermi_Energy, temperature, 1.0, dy, dz, 1, ny, nz)
        {
            this.dE = dE;

            // set band profile with spin-degeneracy
            this.band_gap = Input_Band_Structure.Expand_BandStructure(band_gap, ny);

            this.acceptor_concentration = acceptor_concentration;
            this.donor_concentration = donor_concentration;

            this.dopent_concentration = (SpinResolved_DoubleMatrix)Input_Band_Structure.Expand_BandStructure(donor_concentration - acceptor_concentration, ny);

            // and set relative dopent energies to the conduction band
            this.donor_energy_below_Ec = band_gap / 2.0 - donor_energy;
            this.acceptor_energy_below_Ec = band_gap / 2.0 - acceptor_energy;
        }

        public SpinResolved_DoubleMatrix Get_OneD_Density(DoubleMatrix conduction_band_energy)
        {
            // spin-resolved density
            SpinResolved_DoubleMatrix density = new SpinResolved_DoubleMatrix(ny, nz);

            // Find conduction band density
            SpinResolved_DoubleMatrix conduction_density = Calculate_Conduction_Band_Density(conduction_band_energy);

            // Find valence band density
            SpinResolved_DoubleMatrix valence_density = Calculate_Valence_Band_Density(conduction_band_energy);

            // Calculate donor occupation probability
            SpinResolved_DoubleMatrix acceptor_density = Calculate_Dopent_Density(conduction_band_energy, acceptor_concentration, acceptor_energy_below_Ec);
            SpinResolved_DoubleMatrix donor_density = Calculate_Dopent_Density(conduction_band_energy, donor_concentration, donor_energy_below_Ec);

            // return total density (for 1D; hence the divide-by-lattice-spacing)
            return -1.0 * q_e * (dopent_concentration + valence_density + acceptor_density - conduction_density - donor_density);
        }

        /// <summary>
        /// calculates the conduction band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        SpinResolved_DoubleMatrix Calculate_Conduction_Band_Density(DoubleMatrix conduction_band_energy)
        {
            SpinResolved_DoubleMatrix result = new SpinResolved_DoubleMatrix(ny, nz);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((no_kB_T_above_Ef * kB * temperature - conduction_band_energy.Min()) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    for (int E_step = 0; E_step < no_energy_steps; E_step++)
                    {
                        double energy_above_eF = conduction_band_energy[i, j] + E_step * dE;

                        result[i, j, Spin.Up] += Get_3D_DensityofStates(conduction_band_energy[i, j], energy_above_eF) * Get_Fermi_Function(energy_above_eF);
                        result[i, j, Spin.Down] += Get_3D_DensityofStates(conduction_band_energy[i, j], energy_above_eF) * Get_Fermi_Function(energy_above_eF);
                    }

            return result;
        }

        /// <summary>
        /// calculates the valence band profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        SpinResolved_DoubleMatrix Calculate_Valence_Band_Density(DoubleMatrix conduction_band_energy)
        {
            SpinResolved_DoubleMatrix result = new SpinResolved_DoubleMatrix(ny, nz);
            DoubleMatrix valence_band = conduction_band_energy - band_gap;

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((valence_band.Max() - no_kB_T_above_Ef * kB * temperature) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    for (int E_step = 0; E_step < no_energy_steps; E_step++)
                    {
                        double energy_below_eF = valence_band[i, j] - E_step * dE;

                        // energies are negative (for density of states) as we are considering holes
                        // also, we subtract from "result" as holes are negatively charged
                        result[i, j, Spin.Up] += Get_3D_DensityofStates(-1.0 * valence_band[i, j], -1.0 * energy_below_eF) * (1.0 - Get_Fermi_Function(energy_below_eF));
                        result[i, j, Spin.Down] += Get_3D_DensityofStates(-1.0 * valence_band[i, j], -1.0 * energy_below_eF) * (1.0 - Get_Fermi_Function(energy_below_eF));
                    }

            return result;
        }

        /// <summary>
        /// calculates the dopent density profile (spin-resolved) as a function of depth using a given conduction band potential
        /// </summary>
        SpinResolved_DoubleMatrix Calculate_Dopent_Density(DoubleMatrix conduction_band_energy, DoubleVector dopent_concentration, DoubleVector dopent_energy_below_Ec)
        {
            SpinResolved_DoubleMatrix result = new SpinResolved_DoubleMatrix(ny, nz);

            // calculates the occupation of the dopents 
            // (note that the dopent concentrations are only a function of depth)
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    result[i, j, Spin.Up] = dopent_concentration[j] * Get_Fermi_Function(conduction_band_energy[i, j] - dopent_energy_below_Ec[j]);
                    result[i, j, Spin.Down] = dopent_concentration[j] * Get_Fermi_Function(conduction_band_energy[i, j] - dopent_energy_below_Ec[j]);
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
