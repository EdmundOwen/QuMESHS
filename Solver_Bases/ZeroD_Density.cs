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
    public class ZeroD_Density : Density_Base
    {
        double band_gap;
        double acceptor_conc, donor_conc;
        double acceptor_energy, donor_energy;
        double max_density = 1.0;

        double electron_mass;
        double hole_mass;

        /// <summary>
        /// this value is used as temporary storage of the chemical potential for the root-finding method in order
        /// to integrate up the densities with respect to the energy... it should *ALWAYS* be reset to double.MaxValue
        /// before the method exits!!
        /// </summary>
        double chem_pot = double.MaxValue;

        double no_kB_T = 50.0;

        public ZeroD_Density(ILayer layer, double temperature)
            : base(temperature)
        {
            this.band_gap = layer.Band_Gap;
            this.acceptor_conc = layer.Acceptor_Conc; this.donor_conc = layer.Donor_Conc;
            this.acceptor_energy = layer.Acceptor_Energy; this.donor_energy = layer.Donor_Energy;

            this.electron_mass = layer.Electron_Mass;
            this.hole_mass = layer.Hole_Mass;
        }

        public double Get_ChargeDensity(double mu)
        {
            this.mu = mu;

            // calculate the densities due to the various components for a given chemical potential
            double conductance_electrons = Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), 0.5 * band_gap, Math.Max(0.5 * band_gap, mu) + no_kB_T * Physics_Base.kB * temperature);
            // factor of 2.0 here is for spin degeneracy of the dopents
            // also, exponential factor for donors is (E_d - mu)
            double donor_electrons = donor_conc * 2.0 * Get_Dopent_Fermi_Function(donor_energy, mu, temperature);
            double donor_nuclei = donor_conc;
            double acceptor_nuclei = acceptor_conc;
            // and exponential factor for acceptors is (mu - E_a)
            double acceptor_holes = acceptor_conc * 2.0 * Get_Dopent_Fermi_Function(-acceptor_energy, -mu, temperature);
            double valence_holes = Perform_Integral(new Func<double, double>(Get_Valence_Electron_Integrand), Math.Min(-0.5 * band_gap, mu) - no_kB_T * Physics_Base.kB * temperature, -0.5 * band_gap);

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
            this.mu = mu;

            double donor_density = Get_DonorDensity();
            double acceptor_density = Get_AcceptorDensity();

            return Physics_Base.q_e * (donor_density - acceptor_density);
        }

        double Get_DonorDensity()
        {
            // calculate the densities due to the donors for a given chemical potential
            // factor of 2.0 here is for spin degeneracy of the dopents
            // exponential factor for donors is (E_d - mu)
            double donor_electrons = donor_conc * 2.0 * Get_Dopent_Fermi_Function(donor_energy, mu, temperature);
            double donor_nuclei = donor_conc;

            return donor_nuclei - donor_electrons;
        }

        double Get_AcceptorDensity()
        {
            // calculate the densities due to the acceptors for a given chemical potential
            // factor of 2.0 here is for spin degeneracy of the dopents
            double acceptor_nuclei = acceptor_conc;
            // exponential factor for acceptors is (mu - E_a)
            double acceptor_holes = acceptor_conc * 2.0 * Get_Dopent_Fermi_Function(-acceptor_energy, -mu, temperature);

            return acceptor_nuclei - acceptor_holes;
        }

        public double Get_DopentDensityDeriv(double mu)
        {
            this.mu = mu;

            return -2.0 * Physics_Base.q_e * Physics_Base.q_e * (donor_conc * Get_Dopent_Fermi_Function_Derivative(donor_energy, mu, temperature)
                                                + acceptor_conc * Get_Dopent_Fermi_Function_Derivative(-acceptor_energy, -mu, temperature));
        }

        public double Get_CarrierDensity(double mu)
        {
            this.mu = mu;

            // calculate the densities due to the various components for a given chemical potential
            double conductance_electrons =  Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), 0.5 * band_gap, mu + no_kB_T * Physics_Base.kB * temperature);
    //        if (0.5 * band_gap < mu - no_kB_T * Physics_Base.kB * temperature)
    //        {
    //            double interval = mu - no_kB_T * Physics_Base.kB * temperature - 0.5 *  band_gap;
    //            conductance_electrons += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), mu, 0.5 * band_gap, 0.5 * band_gap + interval);
    //            conductance_electrons += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), mu, 0.5 * band_gap + interval, mu - no_kB_T * Physics_Base.kB * temperature);
    //        }
    //       conductance_electrons += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), mu, mu - no_kB_T * Physics_Base.kB * temperature, mu - 0.5 * no_kB_T * Physics_Base.kB * temperature);
    //       conductance_electrons += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), mu, mu - 0.5 * no_kB_T * Physics_Base.kB * temperature, mu);
    //       conductance_electrons += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), mu, mu, mu + 0.5 * no_kB_T * Physics_Base.kB * temperature);
    //       conductance_electrons += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Integrand), mu, mu + 0.5 * no_kB_T * Physics_Base.kB * temperature, mu + no_kB_T * Physics_Base.kB * temperature);


            // double valence_holes = 0.0;
            //valence_holes += Perform_Integral(new Func<double, double>(Get_Valence_Electron_Integrand), mu, mu - no_kB_T * Physics_Base.kB * temperature, mu);
            //valence_holes += Perform_Integral(new Func<double, double>(Get_Valence_Electron_Integrand), mu, mu, mu + no_kB_T * Physics_Base.kB * temperature);
            double valence_holes = Perform_Integral(new Func<double, double>(Get_Valence_Electron_Integrand), mu - no_kB_T * Physics_Base.kB * temperature, -0.5 * band_gap);

            double density = Physics_Base.q_e *  (valence_holes - conductance_electrons);

            //if (density > max_density)
            //    return max_density;
            //else if (-1.0 * density > max_density)
            //    return -1.0 * max_density;
            //else
                return density;
        }

        public double Get_CarrierDensityDeriv(double mu)
        {
            this.mu = mu;

            double conductance_electrons_deriv = Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), 0.5 * band_gap, mu + no_kB_T * Physics_Base.kB * temperature);// 0.0;
         //   double valence_holes_deriv = 0.0;
            // calculate the densities due to the various components for a given chemical potential
            // only integrate from -nokb kb T to + nokb kb T with the if statement for ensuring no integration in the band gap
            // this is due to the derivative of the fermi function being strongly peaked at the chemical potential
        //    conductance_electrons_deriv += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), mu, mu - no_kB_T * Physics_Base.kB * temperature, mu - 0.5 * no_kB_T * Physics_Base.kB * temperature);
        //    conductance_electrons_deriv += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), mu, mu - 0.5 * no_kB_T * Physics_Base.kB * temperature, mu);
        //    conductance_electrons_deriv += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), mu, mu, mu + 0.5 * no_kB_T * Physics_Base.kB * temperature);
        //    conductance_electrons_deriv += Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), mu, mu + 0.5 * no_kB_T * Physics_Base.kB * temperature, mu + no_kB_T * Physics_Base.kB * temperature);

        //    if (mu - no_kB_T * Physics_Base.kB * temperature < 0.5 * band_gap)
        //        conductance_electrons_deriv = Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), mu, 0.5 * band_gap, Math.Max(0.5 * band_gap, mu) + no_kB_T * Physics_Base.kB * temperature);
        //    else
        //        conductance_electrons_deriv = Perform_Integral(new Func<double, double>(Get_Conduction_Electron_Derivative_Integrand), mu, mu - no_kB_T * Physics_Base.kB * temperature, mu + no_kB_T * Physics_Base.kB * temperature);


            double valence_holes_deriv = Perform_Integral(new Func<double, double>(Get_Valence_Electron_Derivative_Integrand), mu - no_kB_T * Physics_Base.kB * temperature, mu + no_kB_T * Physics_Base.kB * temperature);
            //valence_holes_deriv += Perform_Integral(new Func<double, double>(Get_Valence_Electron_Derivative_Integrand), mu, mu, mu + no_kB_T * Physics_Base.kB * temperature);

            // and for holes
            // only integrate from -nokb kb T to + nokb kb T with the if statement for ensuring no integration in the band gap
      //      if (mu + no_kB_T * Physics_Base.kB * temperature > -0.5 * band_gap)
      //          valence_holes_deriv = Perform_Integral(new Func<double, double>(Get_Valence_Electron_Derivative_Integrand), mu, Math.Min(-0.5 * band_gap, mu) - no_kB_T * Physics_Base.kB * temperature, -0.5 * band_gap);
      //      else
      //          valence_holes_deriv = Perform_Integral(new Func<double, double>(Get_Valence_Electron_Derivative_Integrand), mu, mu - no_kB_T * Physics_Base.kB * temperature, mu + no_kB_T * Physics_Base.kB * temperature);

            return Physics_Base.q_e * Physics_Base.q_e * (valence_holes_deriv - conductance_electrons_deriv);
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
        double Perform_Integral(Func<double, double> function, double lower_limit, double upper_limit)
        {
            if (chem_pot != double.MaxValue)
                throw new Exception("Error - Something isn't reseting the chemical potential... this is very unsafe");
            chem_pot = mu;

            // integrates the density for a given energy up to a given number of kB * T above the conduction band edge
            OneVariableFunction integrand = new OneVariableFunction(function);
            GaussKronrodIntegrator tmp = new GaussKronrodIntegrator();
            double result = tmp.Integrate(integrand, lower_limit, upper_limit);

            // reset chem_pot
            chem_pot = double.MaxValue;

            return result;
        }
        /*
        double test()
        {
            OneVariableFunction integrand = new OneVariableFunction(function);
            GaussKronrodIntegrator gk = new GaussKronrodIntegrator();
            double result = gk.Integrate(integrand, 500.0, 1000.0);

            if (!gk.ToleranceMet)
                throw new Exception();

            return result;
        }

        double function(double x)
        {
            if (x < 712.0)
                return 0.0;

            return Math.Pow((x - 712.0), 0.5) / (Math.Exp((x - 700.0) / 6.0) + 1.0); 
        }
        */
        double Get_Conduction_Electron_Integrand(double energy)
        {
            return Physics_Base.Get_3D_DensityofStates(energy, 0.5 * band_gap, Carrier.electron, electron_mass) * Get_Fermi_Function(energy, chem_pot, temperature);
        }

        double Get_Conduction_Electron_Derivative_Integrand(double energy)
        {
            return Physics_Base.Get_3D_DensityofStates(energy, 0.5 * band_gap, Carrier.electron, electron_mass) * Get_Fermi_Function_Derivative(energy, chem_pot, temperature);
        }

        double Get_Valence_Electron_Integrand(double energy)
        {
            // note that the 3D density of states integration is inverted as the density of states increases with decreasing energy
            return Physics_Base.Get_3D_DensityofStates(energy, -0.5 * band_gap, Carrier.hole, hole_mass) * (1.0 - Get_Fermi_Function(energy, chem_pot, temperature));
        }

        double Get_Valence_Electron_Derivative_Integrand(double energy)
        {
            // note that the 3D density of states integration is inverted as the density of states increases with decreasing energy
            return -1.0 * Physics_Base.Get_3D_DensityofStates(energy, -0.5 * band_gap, Carrier.hole, hole_mass) * Get_Fermi_Function_Derivative(energy, chem_pot, temperature);
        }

        #region unused_methods
        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }
        #endregion
    }
}
