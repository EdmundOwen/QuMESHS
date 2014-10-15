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

        double tx, ty;

        public TwoD_DFTSolver(IExperiment exp)
            : base(exp)
        {
            tx = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
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
            Band_Data band_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref band_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, charge_density, band_pot);
            DoubleHermitianEigDecompServer eig_server = new DoubleHermitianEigDecompServer();
            eig_server.ComputeEigenValueRange(band_pot.mat.Min(), no_kb_T * Physics_Base.kB * temperature);
            eig_server.ComputeVectors = true;
            DoubleHermitianEigDecomp eig_decomp = eig_server.Factor(hamiltonian);

            int max_wavefunction = 0;
            if (eig_decomp.EigenValues.Length != 0)
            {
                double min_eigval = eig_decomp.EigenValues.Min();
                max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;
            }

            DoubleMatrix dens_up = new DoubleMatrix(nx, ny, 0.0);
            DoubleMatrix dens_down = new DoubleMatrix(nx, ny, 0.0);

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double dens_val = 0.0;

                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);
                    }

                    // just share the densities (there is no spin polarisation)
                    dens_up[i, j] = 0.5 * dens_val;
                    dens_down[i, j] = 0.5 * dens_val;
                }

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        /// <summary>
        /// calculates a density of states which assumes that the wave vector is quantised by dk in the transport direction
        /// this gives the solution assuming that the electrons are in a ring
        /// </summary>
        double Get_Discrete_OneD_DoS(double E_c, double no_kb_T)
        {
            double dk = 0.01;
            bool max_j_found = false;
            int j = 0;

            // calculate the number of wave vectors which are below no_kb_T * kB * T above the chemical potential
            while (!max_j_found)
            {
                double E = E_c + Physics_Base.hbar * Physics_Base.hbar * j * j * dk * dk / (2.0 * Physics_Base.mass);
                if (E < no_kb_T * Physics_Base.kB * temperature)
                    j++;
                else
                    max_j_found = true;
            }
            int j_max = j;

            // sum the g(E_j) * f(E_j) contributions for each wave vector k_j = j * dk up to j_max
            double result = 0;
            for (int i = 0; i < j_max; i++)
            {
                double E = Physics_Base.hbar * Physics_Base.hbar * i * i * dk * dk / (2.0 * Physics_Base.mass);
                result += (2.0 * dk / Math.PI) / (Math.Exp((E + E_c) / (Physics_Base.kB * temperature)) + 1.0);
            }

            return result;
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
        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data charge_density_deriv, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data band_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref band_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, charge_density_deriv, band_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            DoubleMatrix dens_up_deriv = new DoubleMatrix(nx, ny, 0.0);
            DoubleMatrix dens_down_deriv = new DoubleMatrix(nx, ny, 0.0);

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double dens_val = 0.0;
                    double dens_val_deriv = 0.0;
                    double result = 0.0;

                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);
                        dens_val_deriv += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS_Deriv(eig_decomp.EigenValue(k), no_kb_T);
                    }

                    result = Physics_Base.q_e * dens_val_deriv;
                   // if (this.alpha_dft != 0.0)
                   //     result += dens_val_deriv * dens_val_deriv * Physics_Base.Get_XC_Potential_Deriv(dens_val);

                    // just share the densities (there is no spin polarisation)
                    dens_up_deriv[i, j] = 0.5 * result;
                    dens_down_deriv[i, j] = 0.5 * result;
                }

                charge_density_deriv = new SpinResolved_Data(new Band_Data(dens_up_deriv), new Band_Data(dens_down_deriv));
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
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(nx * ny);

            // set off diagonal elements 
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // coupling sites in the transverse direction
                    if (i != 0)
                        result[i * ny + j, i * ny + j - ny] = tx; 
                    if (i != nx - 1)
                        result[i * ny + j, i * ny + j + ny] = tx;
                    // coupling sites in the growth direction
                    if (j != 0)
                        result[i * ny + j, i * ny + j - 1] = ty;
                    if (j != ny - 1)
                        result[i * ny + j, i * ny + j + 1] = ty;
                }

            double[,] potential = new double[nx, ny];
            // set diagonal elements
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
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
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(nx * ny);

            // set off diagonal elements 
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // coupling sites in the transverse direction
                    if (i != 0)
                        result[i * ny + j, i * ny + j - ny] = tx; 
                    if (i != nx - 1)
                        result[i * ny + j, i * ny + j + ny] = tx;
                    // coupling sites in the growth direction
                    if (j != 0)
                        result[i * ny + j, i * ny + j - 1] = ty;
                    if (j != ny - 1)
                        result[i * ny + j, i * ny + j + 1] = ty;
                }

            double[,] potential = new double[nx, ny];
            // set diagonal elements
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    potential[i, j] = pot.mat[i, j] + dft_pot.mat[i, j];
                    result[i * ny + j, i * ny + j] = -2.0 * tx + -2.0 * ty + potential[i, j];
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
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    sw.Write(potential[i, j].ToString() + '\t');
                    if (j == ny - 1)
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
                DoubleMatrix tmp = new DoubleMatrix(nx, ny, 0.0);

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        // do not add anything to the density if on the edge of the domain
                        if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                            continue;

                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        tmp[i, j] = -1.0 * Physics_Base.q_e * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);    
                    }

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        sw.Write(tmp[i, j].ToString() + '\t');
                        if (j == ny - 1)
                            sw.WriteLine();
                    }

                sw.Close();
            }
        }
    }
}
