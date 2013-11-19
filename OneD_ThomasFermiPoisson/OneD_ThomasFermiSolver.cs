using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Analysis;
using Solver_Bases;
using Solver_Bases.Layers;

namespace OneD_ThomasFermiPoisson
{
    class OneD_ThomasFermiSolver : Density_Base
    {
        /*
        DoubleVector band_gap;
        SpinResolved_DoubleVector dopent_concentration;
        DoubleVector acceptor_concentration, donor_concentration;
        DoubleVector acceptor_energy_above_Ev, donor_energy_below_Ec;

        double no_kB_T_above_Ef = 10.0;

        // redundent integration parameter
        double dE = 1.0;

        /*public OneD_ThomasFermiSolver(DoubleVector band_gap, DoubleVector acceptor_concentration, DoubleVector donor_concentration, DoubleVector acceptor_energy, DoubleVector donor_energy,
                                        double temperature, double dz, int nz) 
            : base(temperature, 1.0, 1.0, dz, 1, 1, nz, 0.0, 0.0, zmin)
        {
            // set band profile with spin-degeneracy
            this.band_gap = band_gap;

            // set dopent concentrations and copy to a vector of double the length for spin-degeneracy
            this.acceptor_concentration = acceptor_concentration;
            this.donor_concentration = donor_concentration;

            this.dopent_concentration = (SpinResolved_DoubleVector)(donor_concentration - acceptor_concentration);

            // and set relative dopent energies to the conduction band
            this.donor_energy_below_Ec = band_gap / 2.0 - donor_energy;
            this.acceptor_energy_above_Ev = band_gap / 2.0 - acceptor_energy;
        }*/

        public OneD_ThomasFermiSolver(double temperature, double dz, int nz, double zmin)
            : base(temperature, 1.0, 1.0, dz, 1, 1, nz, 0.0, 0.0, zmin)
        {
        }

        /*public SpinResolved_DoubleVector Get_OneD_Density(DoubleVector conduction_band_energy)
        {
            // Find conduction band density
            SpinResolved_DoubleVector conduction_density = Calculate_Conduction_Band_Density(conduction_band_energy);

            // Find valence band density
            DoubleVector valence_band = conduction_band_energy - band_gap;
            SpinResolved_DoubleVector valence_density = Calculate_Valence_Band_Density(valence_band);

            // Calculate donor occupation probability
            SpinResolved_DoubleVector acceptor_density = Calculate_Dopent_Density(- valence_band, acceptor_concentration, acceptor_energy_above_Ev);
            SpinResolved_DoubleVector donor_density = Calculate_Dopent_Density(conduction_band_energy, donor_concentration, donor_energy_below_Ec);

            // return total density (for 1D; hence the divide-by-lattice-spacing)
            return -1.0 * Physics_Base.q_e * (dopent_concentration + valence_density + acceptor_density - conduction_density - donor_density);
        }

        public SpinResolved_DoubleVector Get_OneD_ChargeDensity(DoubleVector chem_pot)
        {
            // spin-resolved charge density
            SpinResolved_DoubleVector charge_density = new SpinResolved_DoubleVector(nz);

            for (int i = 0; i < nz; i++)
            {
                // calculate the charge density at the given point
                ZeroD_Density charge_calc = new ZeroD_Density(band_gap[i], acceptor_concentration[i], acceptor_energy_above_Ev[i], donor_concentration[i], donor_energy_below_Ec[i], temperature);
                double local_chargedensity = charge_calc.Get_ChargeDensity(chem_pot[i]);

                // as there is no spin dependence in this problem yet, just divide the charge into spin-up and spin-down components equally
                charge_density.Spin_Down[i] = 0.5 * local_chargedensity;
                charge_density.Spin_Up[i] = 0.5 * local_chargedensity;
            }

            return charge_density;
        }
        */

        public void Get_OneD_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
                for (int i = 0; i < nz; i++)
                {
                    double z = dz * i + zmin;

                    // get the relevant layer and if it's frozen out, don't recalculate the charge
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z);
                    if (current_Layer.Dopents_Frozen_Out(temperature))
                        continue;

                    // calculate the density at the given point
                    ZeroD_Density charge_calc = new ZeroD_Density(current_Layer, temperature);
                    double local_charge_density = charge_calc.Get_ChargeDensity(chem_pot.vec[i]);

                    // as there is no spin dependence in this problem yet, just divide the charge into spin-up and spin-down components equally
                    density.Spin_Down.vec[i] = 0.5 * local_charge_density;
                    density.Spin_Up.vec[i] = 0.5 * local_charge_density;
                }
        }

        /*
        /// <summary>
        /// calculates the conduction band profile (spin-resolved) as a function of depth using a given conduction band energy
        /// </summary>
        SpinResolved_DoubleVector Calculate_Conduction_Band_Density(DoubleVector conduction_band_energy)
        {
            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(nz);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((no_kB_T_above_Ef * Physics_Base.kB * temperature - conduction_band_energy.Min()) / dE);
            else
                no_energy_steps = 100;

            for (int j = 0; j < nz; j++)
                for (int i = 0; i < no_energy_steps; i++)
                {
                    double energy_above_eF = conduction_band_energy[j] + i * dE;

                    result[j, Spin.Up] += Physics_Base.Get_Electron_3D_DensityofStates(energy_above_eF, conduction_band_energy[j]) * Get_Fermi_Function(energy_above_eF) * dE;
                    result[j, Spin.Down] += Physics_Base.Get_Electron_3D_DensityofStates(energy_above_eF, conduction_band_energy[j]) * Get_Fermi_Function(energy_above_eF) * dE;
                }

            return result;
        }

        /// <summary>
        /// calculates the valence band profile (spin-resolved) as a function of depth using a given conduction band energy
        /// </summary>
        SpinResolved_DoubleVector Calculate_Valence_Band_Density(DoubleVector valence_band)
        {
            SpinResolved_DoubleVector result = new SpinResolved_DoubleVector(nz);

            // Integrate from the minimum point on the conduction band to a given number of kB*T above the fermi surface
            int no_energy_steps;
            if (temperature != 0.0)
                no_energy_steps = (int)((valence_band.Max() + no_kB_T_above_Ef * Physics_Base.kB * temperature) / dE);
            else
                no_energy_steps = 100;

            for (int i = 0; i < no_energy_steps; i++)
                for (int j = 0; j < nz; j++)
                {
                    double energy_below_eF = valence_band[j] - i * dE;

                    // energies are negative (for density of states) as we are considering holes
                    // also, we subtract from "result" as holes are negatively charged
                    result[j, Spin.Up] += Physics_Base.Get_Hole_3D_DensityofStates(-1.0 * energy_below_eF, -1.0 * valence_band[j]) * (1.0 - Get_Fermi_Function(energy_below_eF)) * dE;
                    result[j, Spin.Down] += Physics_Base.Get_Hole_3D_DensityofStates(-1.0 * energy_below_eF, -1.0 * valence_band[j]) * (1.0 - Get_Fermi_Function(energy_below_eF)) * dE;
                }

            return result;
        }

        /// <summary>
        /// calculates the dopent density profile (spin-resolved) as a function of depth using a given conduction band energy
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

        public double Get_Chemical_Potential(int location)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(band_gap[location], acceptor_concentration[location], acceptor_energy_above_Ev[location], donor_concentration[location], donor_energy_below_Ec[location], temperature);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }
        */
          
        public double Get_Chemical_Potential(ILayer[] layers, double z)
        {
            return Get_Chemical_Potential(layers, z, temperature);
        }

        public double Get_Chemical_Potential(ILayer[] layers, double z, double temperature_input)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z), temperature_input);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }
    }
}
