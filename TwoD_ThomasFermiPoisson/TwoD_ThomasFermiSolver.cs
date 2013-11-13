using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_ThomasFermiSolver : Density_Base
    {
        DoubleVector band_gap;
        SpinResolved_DoubleMatrix dopent_concentration;
        DoubleVector acceptor_concentration, donor_concentration;
        DoubleVector acceptor_energy_above_Ev, donor_energy_below_Ec;

        // redundant integration parameters
        //double dE = 1.0;
        //double no_kB_T_above_Ef = 100.0;

        public TwoD_ThomasFermiSolver(DoubleVector band_gap, DoubleVector acceptor_concentration, DoubleVector donor_concentration, DoubleVector acceptor_energy, DoubleVector donor_energy,
                                        double temperature, double dy, double dz, int ny, int nz) 
            : base(temperature, 1.0, dy, dz, 1, ny, nz)
        {
            // set band profile with spin-degeneracy
            this.band_gap = band_gap;// Input_Band_Structure.Expand_BandStructure(band_gap, ny).mat;

            this.acceptor_concentration = acceptor_concentration;
            this.donor_concentration = donor_concentration;

            this.dopent_concentration = (SpinResolved_DoubleMatrix)Input_Band_Structure.Expand_BandStructure(donor_concentration - acceptor_concentration, ny).mat;

            // and set relative dopent energies to the conduction band
            this.donor_energy_below_Ec = band_gap / 2.0 - donor_energy;
            this.acceptor_energy_above_Ev = band_gap / 2.0 - acceptor_energy;
        }

        public SpinResolved_DoubleMatrix Get_TwoD_ChargeDensity(DoubleMatrix chem_pot)
        {
            // spin-resolved density
            SpinResolved_DoubleMatrix charge_density = new SpinResolved_DoubleMatrix(ny, nz);

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    // calculate the density at the given point
                    ZeroD_Density charge_calc = new ZeroD_Density(band_gap[j], acceptor_concentration[j], acceptor_energy_above_Ev[j], donor_concentration[j], donor_energy_below_Ec[j], temperature);
                    double local_charge_density = charge_calc.Get_ChargeDensity(chem_pot[i, j]);

                    // as there is no spin dependence in this problem yet, just divide the charge into spin-up and spin-down components equally
                    charge_density.Spin_Down[i, j] = 0.5 * local_charge_density;
                    charge_density.Spin_Up[i, j] = 0.5 * local_charge_density;
                }

            return charge_density;
        }

        public double Get_Chemical_Potential(int location)
        {
            return Get_Chemical_Potential(location, temperature);
        }

        public double Get_Chemical_Potential(int location, double temperature_input)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(band_gap[location], acceptor_concentration[location], acceptor_energy_above_Ev[location], donor_concentration[location], donor_energy_below_Ec[location], temperature_input);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }

        /*
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
            return -1.0 * Physics_Base.q_e * (dopent_concentration + valence_density + acceptor_density - conduction_density - donor_density);
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
                no_energy_steps = (int)((no_kB_T_above_Ef * Physics_Base.kB * temperature - conduction_band_energy.Min()) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    for (int E_step = 0; E_step < no_energy_steps; E_step++)
                    {
                        double energy_above_eF = conduction_band_energy[i, j] + E_step * dE;

                        result[i, j, Spin.Up] += Physics_Base.Get_Electron_3D_DensityofStates(energy_above_eF, conduction_band_energy[i, j]) * Get_Fermi_Function(energy_above_eF);
                        result[i, j, Spin.Down] += Physics_Base.Get_Electron_3D_DensityofStates(energy_above_eF, conduction_band_energy[i, j]) * Get_Fermi_Function(energy_above_eF);
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
                no_energy_steps = (int)((valence_band.Max() - no_kB_T_above_Ef * Physics_Base.kB * temperature) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    for (int E_step = 0; E_step < no_energy_steps; E_step++)
                    {
                        double energy_below_eF = valence_band[i, j] - E_step * dE;

                        // energies are negative (for density of states) as we are considering holes
                        // also, we subtract from "result" as holes are negatively charged
                        result[i, j, Spin.Up] += Physics_Base.Get_Hole_3D_DensityofStates(-1.0 * energy_below_eF, -1.0 * valence_band[i, j]) * (1.0 - Get_Fermi_Function(energy_below_eF));
                        result[i, j, Spin.Down] += Physics_Base.Get_Hole_3D_DensityofStates(-1.0 * energy_below_eF, -1.0 * valence_band[i, j]) * (1.0 - Get_Fermi_Function(energy_below_eF));
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
        */
    }
}
