/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using Iterative_Greens_Function_Test;

namespace ThreeD_SchrodingerPoissonSolver
{
    class ThreeD_EffectiveBandSolver : ThreeD_Density_Base
    {
        DoubleVector energies;

        double no_kb_T = 50.0;          // number of kb_T to integrate to

        double tx, ty, tz;

        public ThreeD_EffectiveBandSolver(IExperiment exp)
            : base(exp)
        {
            tx = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dx * dx);
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dy * dy);
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dz * dz);
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);
            double[,] xy_energy = new double[nx, ny];

            xy_energy = Solve_Eigenvector_Problem(dft_pot, ref charge_density);

            /*int max_wavefunction = (from val in eig_decomp.EigenValues
                                   where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            double[,] dens_xy = new double[nx, ny];
            // and generate a density for the y-direction
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < max_wavefunction; k++)
                        dens_xy[i, j] += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_Fermi_Function(energies[k]);
*/


            IExperiment exp_green;
            exp_green = new Iterative_Greens_Function_Test.Experiment();
            Iterative_Greens_Function_Test.Iterative_Greens_Function iter = new Iterative_Greens_Function_Test.Iterative_Greens_Function(exp_green, xy_energy);
            double[,] dens_xy = new double[nx, ny];
            double max_energy = Physics_Base.kB * temperature * no_kb_T;
            DoubleMatrix xy_energy_mat = new DoubleMatrix(xy_energy);
            double min_energy = xy_energy_mat.Min();
            double dE = 0.1;
            //double min_energy = -1.0 * max_energy;
            int n_slices = (int)((max_energy - min_energy) / dE);

            for (int n = 1; n < n_slices-1; n++)
            {
                double energy = min_energy + n * dE;
                double[,] dens_per_E = iter.GetDoS(energy);
                double fermi_func = Get_Fermi_Function(energy);

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        dens_xy[i, j] += dens_per_E[i,j] * fermi_func;
                    }
            }

            //double fermi_fact_max = Math.Exp(max_energy / (Physics_Base.kB * temperature)) + 1;
            //double fermi_fact_min = 2.0;
            double[,] dens_per_E_max = iter.GetDoS(max_energy);
            double[,] dens_per_E_min = iter.GetDoS(min_energy);
            double fermi_func_max = Get_Fermi_Function(max_energy);
            double fermi_func_min = Get_Fermi_Function(min_energy);

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    dens_xy[i, j] += 0.5*(dens_per_E_max[i, j] * fermi_func_max + dens_per_E_min[i,j] * fermi_func_min);
                }




            // multiply the z-densities by the xy-density
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                    {
                        charge_density.Spin_Up.vol[k][i, j] *= dens_xy[i, j];
                        charge_density.Spin_Down.vol[k][i, j] *= dens_xy[i, j];
                    }

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = unit_charge * charge_density;
        }

        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data charge_density_deriv, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);
            double[,] xy_energy = new double[nx, ny];

            xy_energy = Solve_Eigenvector_Problem(dft_pot, ref charge_density_deriv);

            /*int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            double[,] dens_xy_deriv = new double[nx, ny];
            // and generate a density for the y-direction
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                for (int k = 0; k < max_wavefunction; k++)
                        dens_xy_deriv[i, j] += DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(k)[i * ny + j]) * Get_Fermi_Function_Derivative(energies[k]);
            */
            IExperiment exp_green;
            exp_green = new Iterative_Greens_Function_Test.Experiment();
            Iterative_Greens_Function_Test.Iterative_Greens_Function iter = new Iterative_Greens_Function_Test.Iterative_Greens_Function(exp_green, xy_energy);
            double[,] dens_xy_deriv = new double[nx, ny];
            double max_energy = Physics_Base.kB * temperature * no_kb_T;
            DoubleMatrix xy_energy_mat = new DoubleMatrix(xy_energy);
            double min_energy = xy_energy_mat.Min();
            double dE = 0.1;
            int n_slices = (int)((max_energy - min_energy) / dE);

            for (int n = 1; n < n_slices - 1; n++)
            {
                double energy = min_energy + n * dE;
                double[,] dens_per_E = iter.GetDoS(energy);
                double fermi_func = Get_Fermi_Function_Derivative(energy);

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        dens_xy_deriv[i, j] += dens_per_E[i, j] * fermi_func;
                    }
            }

            //double fermi_fact_max = Math.Exp(max_energy / (Physics_Base.kB * temperature)) + 1;
            //double fermi_fact_min = 2.0;
            double[,] dens_per_E_max = iter.GetDoS(max_energy);
            double[,] dens_per_E_min = iter.GetDoS(min_energy);
            double fermi_func_max = Get_Fermi_Function_Derivative(max_energy);
            double fermi_func_min = Get_Fermi_Function_Derivative(min_energy);

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    dens_xy_deriv[i, j] += 0.5 * (dens_per_E_max[i, j] * fermi_func_max + dens_per_E_min[i, j] * fermi_func_min);
                }


            // multiply the z-densities by the xy-density
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                    {
                        charge_density_deriv.Spin_Up.vol[k][i, j] *= dens_xy_deriv[i, j];
                        charge_density_deriv.Spin_Down.vol[k][i, j] *= dens_xy_deriv[i, j];
                    }
            
            // and multiply the density derivative by e to get the charge density and by e to convert it to d/dphi (as increasing phi decreases the charge: dn/dphi*-e^2 )
            charge_density_deriv = unit_charge * unit_charge * charge_density_deriv;
        }

        /// <summary>
        /// returns an eigenvector decomposition for the xy-plane after having calculated the density in the
        /// z-direction and storing it in "charge_density"
        /// </summary>
        /// <param name="dft_pot"></param>
        /// <param name="charge_density"></param>
        /// <returns></returns>
        private double[,] Solve_Eigenvector_Problem(Band_Data dft_pot, ref SpinResolved_Data charge_density)
        {
            DoubleHermitianEigDecomp eig_decomp;
            double[,] xy_energy = new double[nx, ny];
            Band_Data dens_up = new Band_Data(nx, ny, nz, 0.0);
            Band_Data dens_down = new Band_Data(nx, ny, nz, 0.0);

            // cycle over xy-plane calculating the ground state energy and eigenstate in the growth direction
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // pull out the chemical potential for this slice
                    double[] z_dft_pot = new double[nz];
                    for (int k = 0; k < nz; k++)
                        z_dft_pot[k] = dft_pot.vol[k][i, j];

                    // calculate its eigendecomposition
                    DoubleHermitianMatrix h_z = Create_Hamiltonian(z_dft_pot, tz, nz);

                    eig_decomp = new DoubleHermitianEigDecomp(h_z);

                    // insert the eigenstate into density and record the local confinement energy
                    xy_energy[i, j] = eig_decomp.EigenValue(0);
                    for (int k = 0; k < nz; k++)
                    {
                        double dens_tot = DoubleComplex.Norm(eig_decomp.EigenVector(0)[k]) * DoubleComplex.Norm(eig_decomp.EigenVector(0)[k]);
                        dens_up.vol[k][i, j] = 0.5 * dens_tot;
                        dens_down.vol[k][i, j] = 0.5 * dens_tot;
                    }
                }

            // calculate the eigenstates in the xy-plane
            /*DoubleHermitianMatrix h_xy = Create_2DEG_Hamiltonian(xy_energy, tx, ty, nx, ny, true, false);
            eig_decomp = new DoubleHermitianEigDecomp(h_xy);
            energies = eig_decomp.EigenValues;*/

            // put the calculated densities into charge_density
            charge_density = new SpinResolved_Data(dens_up, dens_down);

            return xy_energy;
        }

        /// <summary>
        /// returns the eigen-energies for the given potential and charge density
        /// </summary>
        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data pot)
        {
            return energies;
        }

       public Band_Data Get_KS_KE(ILayer[] layers, Band_Data chem_pot)
        {/*
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            // create temporary density for input into...
            Band_Data dens_up = new Band_Data(nx, ny, nz, 0.0);
            Band_Data dens_down = new Band_Data(nx, ny, nz, 0.0);
            SpinResolved_Data dens = new SpinResolved_Data(dens_up, dens_down);

            // calculate eigenvectors
            DoubleHermitianEigDecomp eig_decomp = Solve_Eigenvector_Problem(dft_pot, ref dens);

            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            // generate kinetic energy data
            Band_Data ke = new Band_Data(nx, ny, nz, 0.0);
            Band_Data dens_tot = dens.Spin_Summed_Data;
            double ke_prefactor = -0.5 * Physics_Base.hbar * Physics_Base.hbar / Physics_Base.mass;
            for (int i = 1; i < nx - 1; i++)
                for (int j = 1; j < ny - 1; j++)
                    for (int k = 1; k < nz - 1; k++)
                        for (int n = 0; n < max_wavefunction; n++)
                        {
                            double dy2psi = (dens_tot.vol[k][i, j + 1] + dens_tot.vol[k][i, j - 1] - 2.0 * dens_tot.vol[k][i, j]) * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i]);
                            double dx2psi = DoubleComplex.Norm(eig_decomp.EigenVector(n)[i + 1] + eig_decomp.EigenVector(n)[i - 1] - 2.0 * eig_decomp.EigenVector(n)[i]) * dens_tot.mat[i, j];

                            double psi_div2psi = (dens_tot.mat[i, j] * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i])) * (dx2psi + dy2psi);

                            ke.mat[i, j] += ke_prefactor * psi_div2psi * Get_OneD_DoS(eig_decomp.EigenValue(n), no_kb_T);
                        }
*/
            throw new NotImplementedException();

            //return ke;
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

        DoubleHermitianMatrix Create_2DEG_Hamiltonian(double[,] V, double tx, double ty, int nx, int ny, bool periodic_in_x, bool periodic_in_y)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(nx * ny);

            // set off diagonal elements 
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // coupling sites in the transport direction
                    if (i != 0)
                        result[i * ny + j, i * ny + j - ny] = tx;
                    if (i != nx - 1)
                        result[i * ny + j, i * ny + j + ny] = tx;
                    // coupling sites in the transverse direction
                    if (j != 0)
                        result[i * ny + j, i * ny + j - 1] = ty;
                    if (j != ny - 1)
                        result[i * ny + j, i * ny + j + 1] = ty;
                }

            // set diagonal elements
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    result[i * ny + j, i * ny + j] = -2.0 * tx + -2.0 * ty + V[i, j];

            // and finally, add periodic boundary conditions where appropriate 
            if (periodic_in_x)
                for (int j = 0; j < ny; j++)
                {
                    result[j, j + (nx - 1) * ny] = tx;
                    result[j + (nx - 1) * ny, j] = tx;
                }
            if (periodic_in_y)
                for (int i = 0; i < nx; i++)
                {
                    result[i * ny + 1, i * ny] = ty;
                    result[i * ny, i * ny + 1] = ty;
                }

            return result;
        }

        void Output_Wavefunctions(DoubleHermitianEigDecomp eig, int no_wavefunctions)
        {
            throw new NotImplementedException();
        }
    }
}
