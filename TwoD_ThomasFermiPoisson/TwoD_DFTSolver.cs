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

        public TwoD_DFTSolver(IExperiment exp, Carrier carrier_type)
            : base(exp, carrier_type)
        {
            tx = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dx * dx);
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dy * dy);
        }

        public TwoD_DFTSolver(Experiment exp)
            : this(exp, Carrier.electron)
        {
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

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, dft_pot);
            DoubleHermitianEigDecompServer eig_server = new DoubleHermitianEigDecompServer();
            eig_server.ComputeEigenValueRange(dft_pot.mat.Min(), no_kb_T * Physics_Base.kB * temperature);
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

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double dens_val = 0.0;

                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                    {
                        charge_density.Spin_Up.mat[i, j] = 0.0;
                        charge_density.Spin_Down.mat[i, j] = 0.0;
                        continue;
                    }

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        dens_val += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);
                    }

                    // just share the densities (there is no spin polarisation)
                    charge_density.Spin_Up.mat[i, j] = 0.5 * dens_val;
                    charge_density.Spin_Down.mat[i, j] = 0.5 * dens_val;
                }

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = unit_charge * charge_density;
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
                double E = E_c + Physics_Base.hbar * Physics_Base.hbar * j * j * dk * dk / (2.0 * mass);
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
                double E = Physics_Base.hbar * Physics_Base.hbar * i * i * dk * dk / (2.0 * mass);
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
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, dft_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double dens_val_deriv = 0.0;

                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                    {
                        charge_density_deriv.Spin_Up.mat[i, j] = 0.0;
                        charge_density_deriv.Spin_Down.mat[i, j] = 0.0;
                        continue;
                    }

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 50 * k_b * T above mu = 0
                        dens_val_deriv += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS_Deriv(eig_decomp.EigenValue(k), no_kb_T);
                    }

                    // just share the densities (there is no spin polarisation)
                    charge_density_deriv.Spin_Up.mat[i, j] = 0.5 * dens_val_deriv;
                    charge_density_deriv.Spin_Down.mat[i, j] = 0.5 * dens_val_deriv;
                }

            // and multiply the density derivative by e to get the charge density and by e to convert it to d/dphi (as increasing phi decreases the charge: dn/dphi*-e^2 )
            charge_density_deriv = unit_charge * unit_charge * charge_density_deriv;
        }

        /// <summary>
        /// returns the eigen-energies for the given potential and charge density
        /// </summary>
        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data band_pot = pot.DeepenThisCopy();
            Get_Potential(ref band_pot, layers);

            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, band_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);
            return eig_decomp.EigenValues;
        }

        DoubleHermitianMatrix Create_Hamiltonian(ILayer[] layers, Band_Data pot)
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
                    result[i * ny + j, i * ny + j] = -2.0 * tx + -2.0 * ty + pot.mat[i, j];

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
            Band_Data h_spin_summed = h.Spin_Summed_Data;

            System.IO.StreamWriter sw = new System.IO.StreamWriter(outfile);
            for (int i = 0; i < h_spin_summed.mat.Cols; i++)
                for (int j = 0; j < h_spin_summed.mat.Rows; j++)
                {
                    sw.Write(h_spin_summed.mat[j, i].ToString() + '\t');
                    if (j == h_spin_summed.mat.Rows - 1)
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
                        tmp[i, j] = unit_charge * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_OneD_DoS(eig_decomp.EigenValue(k), no_kb_T);    
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
