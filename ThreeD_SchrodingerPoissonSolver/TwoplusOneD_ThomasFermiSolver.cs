/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

namespace ThreeD_SchrodingerPoissonSolver
{
    class TwoplusOneD_ThomasFermiSolver : ThreeD_Density_Base
    {
        double no_kb_T = 50.0;          // number of kb_T to integrate to

        double tx, ty, tz;
        double g_2d;

        public TwoplusOneD_ThomasFermiSolver(IExperiment exp)
            : base(exp)
        {
            tx = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dx * dx);
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dy * dy);
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dz * dz);

            g_2d = mass / (Math.PI * Physics_Base.hbar);
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            double[,] xy_energy = Solve_Eigenvector_Problem(dft_pot, ref charge_density);
            double beta = 1.0 / (Physics_Base.kB * temperature);

            // multiply the z-densities by the xy-density
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                    {
                        charge_density.Spin_Up.vol[k][i, j] *= g_2d * Math.Log(1.0 + Math.Exp(-1.0 * beta * xy_energy[i, j])) / beta;
                        charge_density.Spin_Down.vol[k][i, j] *= g_2d * Math.Log(1.0 + Math.Exp(-1.0 * beta * xy_energy[i, j])) / beta;
                    }

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = unit_charge * charge_density;
        }

        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data charge_density_deriv, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            double[,] xy_energy = Solve_Eigenvector_Problem(dft_pot, ref charge_density_deriv);
            double beta = 1.0 / (Physics_Base.kB * temperature);

            // multiply the z-densities by the xy-density
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < nz; k++)
                    {
                        charge_density_deriv.Spin_Up.vol[k][i, j] *= g_2d * Get_Fermi_Function(xy_energy[i, j]);
                        charge_density_deriv.Spin_Down.vol[k][i, j] *= g_2d * Get_Fermi_Function(xy_energy[i, j]);
                    }
            
            // and multiply the density derivative by e to get the charge density and by e to convert it to d/dphi (as increasing phi decreases the charge: dn/dphi*-e^2 )
            charge_density_deriv = unit_charge * unit_charge * charge_density_deriv;
        }

        /// <summary>
        /// returns an eigenvector decomposition for the xy-plane after having calculated the density in the
        /// z-direction and storing it in "charge_density"
        /// </summary>
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

            // put the calculated densities into charge_density
            charge_density = new SpinResolved_Data(dens_up, dens_down);

            return xy_energy;
        }

        /// <summary>
        /// returns the eigen-energies for the given potential and charge density
        /// </summary>
        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data pot)
        {
            throw new NotImplementedException();
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
    }
}