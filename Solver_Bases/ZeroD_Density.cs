using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Analysis;
using Solver_Bases.Layers;

namespace Solver_Bases
{
    /// <summary>
    /// a class containing everything necessary to solve the zero dimensional density problem
    /// </summary>
    public class ZeroD_Density
    {
        double band_gap;
        double acceptor_conc, donor_conc;
        double acceptor_energy, donor_energy;
        double temperature;
        double max_density = 1.0;

        /// <summary>
        /// this value is used as temporary storage of the chemical potential for the root-finding method in order
        /// to integrate up the densities with respect to the energy... it should *ALWAYS* be reset to double.MaxValue
        /// before the method exits!!
        /// </summary>
        double chem_pot = double.MaxValue;

        int no_kB_T = 10;

        public ZeroD_Density(double band_gap, double acceptor_conc, double acceptor_energy, double donor_conc, double donor_energy, double temperature)
        {
            // set band profile with spin-degeneracy
            this.band_gap = band_gap;

            // set dopent concentrations and copy to a vector of double the length for spin-degeneracy
            this.acceptor_conc = acceptor_conc;
            this.donor_conc = donor_conc;

            // and set relative dopent energies to the conduction band
            this.donor_energy = donor_energy;
            this.acceptor_energy = acceptor_energy;

            this.temperature = temperature;
        }

        public ZeroD_Density(ILayer layer, double temperature)
        {
            this.band_gap = layer.Band_Gap;
            this.acceptor_conc = layer.Acceptor_Conc; this.donor_conc = layer.Donor_Conc;
            this.acceptor_energy = layer.Acceptor_Energy; this.donor_energy = layer.Donor_Energy;

            this.temperature = temperature;
        }

        public double Get_ChargeDensity(double mu)
        {
            // calculate the densities due to the various components for a given chemical potential
            double conductance_electrons = Get_Conductance_Electron_Density(mu);
            // factor of 2.0 here is for spin degeneracy of the dopents
            // also, exponential factor for donors is (E_d - mu)
            double donor_electrons = donor_conc * 2.0 * Physics_Base.Get_Dopent_Fermi_Function(donor_energy, mu, temperature);//(0.5 * band_gap - donor_energy, mu, temperature);
            double donor_nuclei = donor_conc;
            double acceptor_nuclei = acceptor_conc;
            // and exponential factor for acceptors is (mu - E_a)
            double acceptor_holes = acceptor_conc * 2.0 * Physics_Base.Get_Dopent_Fermi_Function(-acceptor_energy, -mu, temperature);//(-1.0 * (-0.5 * band_gap + acceptor_energy), -mu, temperature);
            double valence_holes = Get_Valence_Electron_Density(mu);

            double density = Physics_Base.q_e * (donor_nuclei - acceptor_nuclei + acceptor_holes + valence_holes - conductance_electrons - donor_electrons);

            if (density > max_density)
                return max_density;
            else if (-1.0 * density > max_density)
                return -1.0 * max_density;
            else
                return density;
        }

        public double Get_DopentDensity(double mu)
        {
            // calculate the densities due to the various components for a given chemical potential
            // factor of 2.0 here is for spin degeneracy of the dopents
            // also, exponential factor for donors is (E_d - mu)
            double donor_electrons = donor_conc * 2.0 * Physics_Base.Get_Dopent_Fermi_Function(donor_energy, mu, temperature);//(0.5 * band_gap - donor_energy, mu, temperature);
            double donor_nuclei = donor_conc;
            double acceptor_nuclei = acceptor_conc;
            // and exponential factor for acceptors is (mu - E_a)
            double acceptor_holes = acceptor_conc * 2.0 * Physics_Base.Get_Dopent_Fermi_Function(-acceptor_energy, -mu, temperature);//(-1.0 * (-0.5 * band_gap + acceptor_energy), -mu, temperature);

            return Physics_Base.q_e * (donor_nuclei - acceptor_nuclei + acceptor_holes - donor_electrons);
        }

        public double Get_CarrierDensity(double mu)
        {
            // calculate the densities due to the various components for a given chemical potential
            double conductance_electrons = Get_Conductance_Electron_Density(mu);
            double valence_holes = Get_Valence_Electron_Density(mu);

            double density = Physics_Base.q_e * (valence_holes - conductance_electrons);

            if (density > max_density)
                return max_density;
            else if (-1.0 * density > max_density)
                return -1.0 * max_density;
            else
                return density;
        }

        public double Get_Equilibrium_Chemical_Potential()
        {
            OneVariableFunction charge_func = new OneVariableFunction(new Func<double, double>(Get_ChargeDensity));

            // initialise root finder and search for a root for the chemical potential
            RiddersRootFinder finder = new RiddersRootFinder();
            double result = finder.Find(charge_func, -1.0 * band_gap, band_gap);

            return result;
        }

        /// <summary>
        /// calculates the density of electrons in the conduction band for a given chemical potential
        /// </summary>
        double Get_Conductance_Electron_Density(double mu)
        {
            if (chem_pot != double.MaxValue)
                throw new Exception("Error - Something isn't reseting the chemical potential... this is very unsafe");
            chem_pot = mu;

            // integrates the density for a given energy up to a given number of kB * T above the conduction band edge
            OneVariableFunction cond_elec_integrand = new OneVariableFunction(new Func<double, double>(Get_Conduction_Electron_Integrand));
            double result = cond_elec_integrand.Integrate(0.5 * band_gap, Math.Max(0.5 * band_gap, mu) + no_kB_T * Physics_Base.kB * temperature);

            // reset chem_pot
            chem_pot = double.MaxValue;

            return result;
        }

        double Get_Conduction_Electron_Integrand(double energy)
        {
            return Physics_Base.Get_Electron_3D_DensityofStates(energy, 0.5 * band_gap) * Physics_Base.Get_Fermi_Function(energy, chem_pot, temperature);
        }

        /// <summary>
        /// calculates the density of holes in the valence band for a given chemical potential
        /// </summary>
        double Get_Valence_Electron_Density(double mu)
        {
            if (chem_pot != double.MaxValue)
                throw new Exception("Error - Something isn't reseting the chemical potential... this is very unsafe");
            chem_pot = mu;

            // integrates the density for a given energy up to a given number of kB * T below the valence band edge
            OneVariableFunction val_elec_integrand = new OneVariableFunction(new Func<double, double>(Get_Valence_Electron_Integrand));
            double result = val_elec_integrand.Integrate(Math.Min(-0.5 * band_gap, mu) - no_kB_T * Physics_Base.kB * temperature, -0.5 * band_gap);

            // reset chem_pot
            chem_pot = double.MaxValue;

            return result;
        }

        double Get_Valence_Electron_Integrand(double energy)
        {
            // note that the 3D density of states integration is inverted as the density of states increases with decreasing energy
            return Physics_Base.Get_Hole_3D_DensityofStates(energy, -0.5 * band_gap) * (1.0 - Physics_Base.Get_Fermi_Function(energy, chem_pot, temperature));
        }
    }
}
