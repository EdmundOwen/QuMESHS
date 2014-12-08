using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using CenterSpace.NMath.Analysis;

namespace OneD_ThomasFermiPoisson
{
    class OneD_DFTSolver : OneD_Density_Base
    {
        double no_kb_T = 50;    // number of kb_T to integrate to
        double t;
        int max_wavefunction = 0;

        public OneD_DFTSolver(IExperiment exp)
            : base(exp)
        {
            t = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // interpolate the input charge density and chemical potential onto a reduced domain to simplify DFT solve
            SpinResolved_Data dft_dens = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));
            Band_Data dft_pot = new Band_Data(new DoubleVector(nz));
            Interpolate_DFT_Grid(ref dft_dens, ref dft_pot, charge_density, chem_pot);
            Get_Potential(ref dft_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, dft_dens, dft_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            DoubleVector dens_up = new DoubleVector(nz, 0.0);
            DoubleVector dens_down = new DoubleVector(nz, 0.0);

            for (int j = 0; j < nz; j++)
            {
                double dens_val = 0.0;
                for (int i = 0; i < max_wavefunction; i++)
                {
                    // and integrate the density of states at this position for this eigenvector from the minimum energy to
                    // (by default) 50 * k_b * T above mu = 0
                    //dens_val += dens_of_states.Integrate(min_eigval, no_kb_T * Physics_Base.kB * temperature);
                    dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(i)[j]) * DoubleComplex.Norm(eig_decomp.EigenVector(i)[j]) * Get_TwoD_DoS(eig_decomp.EigenValue(i));
                }

                // just share the densities (there is no spin polarisation)
                dens_up[j] = 0.5 * dens_val;
                dens_down[j] = 0.5 * dens_val;
            }

            // and multiply the density by -e to get the charge density (as these are electrons)
            dft_dens = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));

            Insert_DFT_Charge(ref charge_density, dft_dens);
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            // artificially deepen the copies of spin up and spin down
            Band_Data tmp_spinup = new Band_Data(new DoubleVector(carrier_density_deriv.Spin_Up.vec.Length));
            Band_Data tmp_spindown = new Band_Data(new DoubleVector(carrier_density_deriv.Spin_Up.vec.Length));

            for (int i = 0; i < carrier_density_deriv.Spin_Up.vec.Length; i++)

                {
                    tmp_spinup.vec[i] = carrier_density_deriv.Spin_Up.vec[i];
                    tmp_spindown.vec[i] = carrier_density_deriv.Spin_Down.vec[i];
                }

            SpinResolved_Data new_density = new SpinResolved_Data(tmp_spinup, tmp_spindown);

            // finally, get the charge density and send it to this new array
            Get_ChargeDensity_Deriv(layers, ref new_density, chem_pot);

            return new_density;
        }

        /// <summary>
        /// Calculate the charge density for this potential
        /// NOTE!!! that the boundary potential is (essentially) set to infty by ringing the density with a set of zeros.
        ///         this prevents potential solvers from extrapolating any residual density at the edge of the eigenstate
        ///         solution out of the charge density calculation domain
        /// </summary>
        /// <param name="layers"></param>
        /// <param name="charge_density_deriv"></param>
        /// <param name="chem_pot"></param>
        public void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data charge_density_deriv, Band_Data chem_pot)
        {
            // interpolate the input charge density and chemical potential onto a reduced domain to simplify DFT solve
            SpinResolved_Data dft_dens = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));
            Get_ChargeDensity(layers, ref charge_density_deriv, chem_pot);
            Band_Data dft_pot = new Band_Data(new DoubleVector(nz));
            Interpolate_DFT_Grid(ref dft_dens, ref dft_pot, charge_density_deriv, chem_pot);
            Get_Potential(ref dft_pot, layers);

       //     // set dft_dens to zero so that the hamiltonian doesn't include the XC term
       //     dft_dens = 0.0 * dft_dens;

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, dft_dens, dft_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            max_wavefunction = (from val in eig_decomp.EigenValues
                                where val < no_kb_T * Physics_Base.kB * temperature
                                select val).ToArray().Length;

            DoubleVector dens_up_deriv = new DoubleVector(nz, 0.0);
            DoubleVector dens_down_deriv = new DoubleVector(nz, 0.0);

                for (int j = 0; j < nz; j++)
                {
                double dens_val = 0.0;
                for (int i = 0; i < max_wavefunction; i++)
                {
                    // and integrate the density of states at this position for this eigenvector from the minimum energy to
                    // (by default) 50 * k_b * T above mu = 0
                    //dens_val += dens_of_states.Integrate(min_eigval, no_kb_T * Physics_Base.kB * temperature);
                    dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(i)[j]) * DoubleComplex.Norm(eig_decomp.EigenVector(i)[j]) * Physics_Base.mass / (2.0 * Math.PI * Physics_Base.hbar * Physics_Base.hbar);// *Physics_Base.Get_Fermi_Function_Derivative(eig_decomp.EigenValue(i), 0.0, no_kb_T) / (Math.PI * Physics_Base.hbar * Physics_Base.hbar);
                }

                // just share the densities (there is no spin polarisation)
                dens_up_deriv[j] = 0.5 * dens_val;
                dens_down_deriv[j] = 0.5 * dens_val;
            }

            // and multiply the density by -e to get the charge density (as these are electrons)
            dft_dens = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up_deriv), new Band_Data(dens_down_deriv));

            Insert_DFT_Charge(ref charge_density_deriv, dft_dens);
        }

        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public DoubleVector Get_EnergyLevels(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // interpolate the input charge density and chemical potential onto a reduced domain to simplify DFT solve
            SpinResolved_Data dft_dens = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));
            Band_Data dft_pot = new Band_Data(new DoubleVector(nz));
            Interpolate_DFT_Grid(ref dft_dens, ref dft_pot, charge_density, chem_pot);
            Get_Potential(ref dft_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, dft_dens, dft_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            System.IO.StreamWriter sw_pot = new System.IO.StreamWriter("tmp_pot.dat");
            System.IO.StreamWriter sw_dens = new System.IO.StreamWriter("tmp_dens.dat");

            Band_Data dft_dens_spin_summed = dft_dens.Spin_Summed_Data;
            for (int i = 0; i < dft_pot.Length; i++)
            {
                sw_dens.WriteLine(dft_dens_spin_summed[i].ToString());
                sw_pot.WriteLine(dft_pot[i].ToString());
            }

            sw_dens.Close();
            sw_pot.Close();

            return eig_decomp.EigenValues;
        }

        /// <summary>
        /// protects carrier_density such that its values are not overriden
        /// </summary>
        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            SpinResolved_Data tmp_dens = carrier_density.DeepenThisCopy();
            Get_ChargeDensity(layers, ref tmp_dens, ref dopent_density, chem_pot);
            return tmp_dens + dopent_density;
        }


        void Insert_DFT_Charge(ref SpinResolved_Data charge_density, SpinResolved_Data dft_dens)
        {
            for (int i = 0; i < nz; i++)
            {
                int init_index = (int)Math.Floor((i * dz - zmin_pot + zmin) / dz_pot); // index for the initial domain (ie. for charge_density and band_offset)

                // no interpolation (it doesn't work...)
                charge_density.Spin_Up[init_index] = dft_dens.Spin_Up[i];
                charge_density.Spin_Down[init_index] = dft_dens.Spin_Down[i];
            }
        }

        void Interpolate_DFT_Grid(ref SpinResolved_Data dft_dens, ref Band_Data dft_band_offset, SpinResolved_Data charge_density, Band_Data band_offset)
        {
            for (int i = 0; i < nz; i++)
            {
                int init_index = (int)Math.Floor((i * dz - zmin_pot + zmin) / dz_pot); // index for the initial domain (ie. for charge_density and band_offset)

                // no interpolation (it doesn't work...)
                dft_dens.Spin_Up[i] = charge_density.Spin_Up[init_index];
                dft_dens.Spin_Down[i] = charge_density.Spin_Down[init_index];
                dft_band_offset[i] = band_offset[init_index];
            }
        }

        void Get_Potential(ref Band_Data dft_band_offset, ILayer[] layers)
        {
            for (int i = 0; i < nz; i++)
            {
                double pos = zmin + i * dz;
                double band_gap = Geom_Tool.GetLayer(layers, pos).Band_Gap;
                dft_band_offset[i] = 0.5 * band_gap - dft_band_offset[i];
            }
        }

        double Get_TwoD_DoS(double tmp_eigval)
        {
            // calculate the density of states integral directly
            double alpha = Physics_Base.mass / (Physics_Base.hbar * Physics_Base.hbar * 2.0 * Math.PI);
            double beta = 1.0 / (Physics_Base.kB * temperature);
            OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double E) => 1.0 / (Math.Exp(beta * E) + 1)));
            dos_integrand.Integrator = new GaussKronrodIntegrator();
            return alpha * dos_integrand.Integrate(tmp_eigval, no_kb_T * Physics_Base.kB * temperature);
        }

        DoubleHermitianMatrix Create_Hamiltonian(ILayer[] layers, SpinResolved_Data charge_density, Band_Data pot)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(nz);

            // set off diagonal elements
            for (int i = 0; i < nz - 1; i++)
            {
                result[i + 1, i] = t; result[i, i + 1] = t;
            }

            double[] potential = new double[nz];
            Band_Data charge_dens_spin_summed = charge_density.Spin_Summed_Data;
            // set diagonal elements
            for (int i = 0; i < nz; i++)
            {
                potential[i] = pot.vec[i] + Physics_Base.Get_XC_Potential(charge_dens_spin_summed.vec[i]);
                result[i, i] = -2.0 * t + potential[i];
            }

            return result;
        }

        public int No_Wavefunctions
        {
            get { return max_wavefunction; }
        }

        double zmin_pot;
        public double Zmin_Pot
        {
            set { zmin_pot = value; }
        }
        double dz_pot;
        public double Dz_Pot
        {
            set { dz_pot = value; }
        }
    }
}
