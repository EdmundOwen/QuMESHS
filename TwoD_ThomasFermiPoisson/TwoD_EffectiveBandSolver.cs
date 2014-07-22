using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;

namespace TwoD_ThomasFermiPoisson
{
    public class TwoD_EffectiveBandSolver : TwoD_Density_Base
    {
        DoubleVector energies;

        double no_kb_T = 50.0;          // number of kb_T to integrate to

        double ty, tz;

        public TwoD_EffectiveBandSolver(Experiment exp)
            : base(exp)
        {
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            DoubleHermitianEigDecomp eig_decomp;

            double[] y_energy = new double[ny];
            DoubleMatrix dens_up = new DoubleMatrix(ny, nz, 0.0);
            DoubleMatrix dens_down = new DoubleMatrix(ny, nz, 0.0);

            // cycle over y-direction calculating the ground state energy and eigenstate in the z-direction
            for (int i = 0; i < ny; i++)
            {
                // pull out the chemical potential for this slice
                double[] z_dft_pot = new double[nz];
                for (int j = 0; j < nz; j++)
                    z_dft_pot[j] = dft_pot.mat[i, j];

                // calculate its eigendecomposition
                DoubleHermitianMatrix h_z = Create_Hamiltonian(z_dft_pot, tz, nz);
                eig_decomp = new DoubleHermitianEigDecomp(h_z);

                // insert the eigenstate into density and record the local confinement energy
                y_energy[i] = eig_decomp.EigenValue(0);
                for (int j = 0; j < nz; j++)
                {
                    dens_up[i, j] = DoubleComplex.Norm(0.5 * eig_decomp.EigenVector(0)[j]) * DoubleComplex.Norm(0.5 * eig_decomp.EigenVector(0)[j]);
                    dens_down[i, j] = DoubleComplex.Norm(0.5 * eig_decomp.EigenVector(0)[j]) * DoubleComplex.Norm(0.5 * eig_decomp.EigenVector(0)[j]);
                }
            }

            // calculate the eigenstates in the y-direction
            DoubleHermitianMatrix h_y = Create_Hamiltonian(y_energy, ty, ny);
            eig_decomp = new DoubleHermitianEigDecomp(h_y);
            energies = eig_decomp.EigenValues;

            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;
            double[] dens_y = new double[ny];

            // and generate a density for the y-direction
            for (int i = 0; i < ny; i++)
                for (int k = 0; k < max_wavefunction; k++)
                    dens_y[i] += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);

            // multiply the z-densities by the y-density
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    dens_up[i, j] *= dens_y[i];
                    dens_down[i, j] *= dens_y[i];
                }
            
            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        /// <summary>
        /// returns the eigen-energies for the given potential and charge density
        /// </summary>
        public DoubleVector Get_EnergyLevels(ILayer[] layers, SpinResolved_Data charge_density, Band_Data pot)
        {
            return energies;
        }
        
        DoubleHermitianMatrix Create_Hamiltonian(double[] V, double t, int N)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(N);

            // set off diagonal elements 
            for (int i = 0; i < N; i++)
                {
                    // coupling sites in the growth direction
                    if (i != 0)
                        result[i - 1, i] = t;
                    if (i != N - 1)
                        result[i + 1, i] = t;
                }

            // set diagonal elements
            for (int i = 0; i < N; i++)
                result[i, i] = -2.0 * t + V[i];
                
            return result;
        }

        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }
    }
}
