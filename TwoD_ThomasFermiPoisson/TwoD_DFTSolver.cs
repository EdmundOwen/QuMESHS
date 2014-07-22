using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;

namespace TwoD_ThomasFermiPoisson
{
    public class TwoD_DFTSolver : TwoD_Density_Base
    {
        double no_kb_T = 50.0;          // number of kb_T to integrate to

        double ty, tz;

        public TwoD_DFTSolver(Experiment exp)
            : base(exp)
        {
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);
        }

        /// <summary>
        /// Calculate the charge density for this potential
        /// NOTE!!! that the boundary potential is (essentially) set to infty by ringing the density with a set of zeros.
        ///         this prevents potential solvers from extrapolating any residual density at the edge of the eigenstate
        ///         solution out of the charge density calculation domain
        /// </summary>
        /// <param name="layers"></param>
        /// <param name="charge_density"></param>
        /// <param name="chem_pot"></param>
        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, charge_density, dft_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            DoubleMatrix dens_up = new DoubleMatrix(ny, nz, 0.0);
            DoubleMatrix dens_down = new DoubleMatrix(ny, nz, 0.0);

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    double dens_val = 0.0;

                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                        continue;

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * nz + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * nz + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);
                    }

                    // just share the densities (there is no spin polarisation)
                    dens_up[i, j] = 0.5 * dens_val;
                    dens_down[i, j] = 0.5 * dens_val;
                }

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        /// <summary>
        /// Calculate the charge density for this potential
        /// NOTE!!! that the boundary potential is (essentially) set to infty by ringing the density with a set of zeros.
        ///         this prevents potential solvers from extrapolating any residual density at the edge of the eigenstate
        ///         solution out of the charge density calculation domain
        /// </summary>
        /// <param name="layers"></param>
        /// <param name="charge_density"></param>
        /// <param name="chem_pot"></param>
        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, charge_density, dft_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            DoubleMatrix dens_up = new DoubleMatrix(ny, nz, 0.0);
            DoubleMatrix dens_down = new DoubleMatrix(ny, nz, 0.0);

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    double dens_val = 0.0;

                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                        continue;

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * nz + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * nz + j]) * Get_OneD_DoS_Deriv(eig_decomp.EigenValue(k), no_kb_T);
                    }

                    // just share the densities (there is no spin polarisation)
                    dens_up[i, j] = 0.5 * dens_val;
                    dens_down[i, j] = 0.5 * dens_val;
                }


            /////////// THERE SHOULD BE A MINUS SIGN HERE!!! NO IDEA WHERE ITS'S GONE /////////////////

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density =  Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        /// <summary>
        /// returns the eigen-energies for the given potential and charge density
        /// </summary>
        public DoubleVector Get_EnergyLevels(ILayer[] layers, SpinResolved_Data charge_density, Band_Data pot)
        {
            Get_Potential(ref pot, layers);
            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, charge_density, pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);
            return eig_decomp.EigenValues;
        }

        DoubleHermitianMatrix Create_Hamiltonian(ILayer[] layers, SpinResolved_Data charge_density, Band_Data pot)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(ny * nz);

            // set off diagonal elements 
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    // coupling sites in the growth direction
                    if (i != 0)
                        result[i * nz + j, i * nz + j - nz] = ty; 
                    if (i != ny - 1)
                        result[i * nz + j, i * nz + j + nz] = ty;
                    // coupling sites in the transverse direction
                    if (j != 0)
                        result[i * nz + j, i * nz + j - 1] = tz;
                    if (j != nz - 1)
                        result[i * nz + j, i * nz + j + 1] = tz;
                }

            double[,] potential = new double[ny, nz];
            // set diagonal elements
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    potential[i, j] = pot.mat[i, j];// +Physics_Base.Get_XC_Potential(charge_density.Spin_Summed_Data.mat[i, j]);  This should already be included in the input chemical potential
                    result[i * nz + j, i * nz + j] = -2.0 * ty + -2.0 * tz + potential[i, j];
                }

            return result;
        }

        public void Write_Out_Hamiltonian(DoubleHermitianMatrix h)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter("hamiltonian.dat");
            for (int i = 0; i < h.Cols; i++)
                for (int j = 0; j < h.Rows; j++)
                {
                    sw.Write(h[i, j].Real.ToString() + '\t');
                    if (j == h.Rows - 1)
                        sw.WriteLine();
                }

            sw.Close();
        }

        public void Write_Out_Density(SpinResolved_Data h, string outfile)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter(outfile);
            for (int i = 0; i < h.Spin_Summed_Data.mat.Cols; i++)
                for (int j = 0; j < h.Spin_Summed_Data.mat.Rows; j++)
                {
                    sw.Write(h.Spin_Summed_Data.mat[j, i].ToString() + '\t');
                    if (j == h.Spin_Summed_Data.mat.Rows - 1)
                        sw.WriteLine();
                }

            sw.Close();
        }

        public void Write_Out_Potential(double[,] potential, string outfile)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter(outfile);
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    sw.Write(potential[i, j].ToString() + '\t');
                    if (j == nz - 1)
                        sw.WriteLine();
                }

            sw.Close();
        }

        void Write_Out_Wavefunction_Densities(DoubleHermitianEigDecomp eig_decomp)
        {
            double min_eigval = eig_decomp.EigenValues.Min();
            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            for (int k = 0; k < max_wavefunction; k++)
            {
                System.IO.StreamWriter sw = new System.IO.StreamWriter("dens_wavefunction_" + k.ToString("00"));
                DoubleMatrix tmp = new DoubleMatrix(ny, nz, 0.0);

                for (int i = 0; i < ny; i++)
                    for (int j = 0; j < nz; j++)
                    {
                        // do not add anything to the density if on the edge of the domain
                        if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                            continue;

                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        tmp[i, j] = -1.0 * Physics_Base.q_e * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * nz + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * nz + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);    
                    }

                for (int i = 0; i < ny; i++)
                    for (int j = 0; j < nz; j++)
                    {
                        sw.Write(tmp[i, j].ToString() + '\t');
                        if (j == nz - 1)
                            sw.WriteLine();
                    }

                sw.Close();
            }
        }
    }
}
