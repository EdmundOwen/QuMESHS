using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using System.IO;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_SO_DFTSolver : TwoD_Density_Base
    {
        double no_kb_T = 50.0;          // number of kb_T to integrate to
        double delta_k = 0.01;           // integration variable in k for band structure
        double eps = 1.0e-9;            // minimum density added before convergence

        DoubleComplex I = new DoubleComplex(0.0, 1.0);

        double tx, ty;
        double theta_x, theta_y;
        double alpha;
        double g_1D;                    // one dimensional density of states in k; g(k) = dn/dk
                                        // NOTE: this is without spin and k degeneracy and is, therefore, a factor of 1/4 times the non-SOI case
        double E_min = -100.0;          // minimum energy value for Lanczos matrix diagonalisation routine in meV
                                        // should be less than minimum energy expected for any k-vector

        Band_Data dV_x, dV_y, dV_xy;
        double kmax = 0.0;

        public TwoD_SO_DFTSolver(Experiment exp)
            : base(exp)
        {
            tx = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            double r_so = 117.1 * (0.01 / Physics_Base.q_e);                                                    // r^{6c6c}_{41} for InAs as reported by Winkler (Table 6.6, p.87)
                                                                                                                // NOTE: the result there is in e A^2... to convert to nm^2, we divide by shown factors
            alpha = r_so / Physics_Base.q_e;
            
            theta_x = r_so * Physics_Base.mass * dx / (Physics_Base.q_e / Physics_Base.hbar);
            theta_y = -1.0 * r_so * Physics_Base.mass * dy / (Physics_Base.q_e / Physics_Base.hbar);

            g_1D = 0.5 / Math.PI;
        }

        public void Get_SOI_parameters(Band_Data chem_pot)
        {
            // initialise matrices with same dimension as the chem_pot data which we will receive later
            dV_x = new Band_Data(chem_pot.mat.Rows, chem_pot.mat.Cols, 0.0);
            dV_y = new Band_Data(chem_pot.mat.Rows, chem_pot.mat.Cols, 0.0);
            dV_xy = new Band_Data(chem_pot.mat.Rows, chem_pot.mat.Cols, 0.0);

            // calculate the derivatives of chem_pot in y and z and the second mixed derivative
            for (int i = 1; i < nx - 1; i++)
                for (int j = 1; j < ny - 1; j++)
                {
                    dV_x.mat[i, j] = (chem_pot.mat[i + 1, j] - chem_pot.mat[i - 1, j]) / (2.0 * dx);
                    dV_y.mat[i, j] = (chem_pot.mat[i, j + 1] - chem_pot.mat[i, j - 1]) / (2.0 * dy);
                    dV_xy.mat[i, j] = (chem_pot.mat[i + 1, j + 1] - chem_pot.mat[i - 1, j + 1] - chem_pot.mat[i + 1, j - 1] + chem_pot.mat[i - 1, j - 1]) / (4.0 * dx * dy);
                }
        }

        /// <summary>
        /// </summary>
        /// <param name="layers"></param>
        /// <param name="charge_density"></param>
        /// <param name="chem_pot"></param>
        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            Get_SOI_parameters(chem_pot);

         //   if (dV_x == null)
         //       throw new Exception("Error - Band structure derivatives are null!  Have you initiated this type properly by calling Get_SOI_parameters(Band_Data chem_pot)?");

            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);
            E_min = dft_pot.Min();

            // reset charge density
            charge_density = 0.0 * charge_density;
            
            // integrate up density from k = 0 to k = k_F
            bool integrated = false;
            int count = 0;
            double k = 0;

            // create initial hamiltonian and eigendecomposition
            int max_wavefunction_old_p, max_wavefunction_old_m;
            DoubleHermitianMatrix hamiltonian_p = Create_H_Using_SOIP(dft_pot, k);
            DoubleHermitianEigDecomp eig_decomp_old_p = Diagonalise_Hamiltonian(hamiltonian_p, out max_wavefunction_old_p);
            DoubleHermitianEigDecomp eig_decomp_old_m = Diagonalise_Hamiltonian(hamiltonian_p, out max_wavefunction_old_m);
            
            while (!integrated)
            {
                // increment the wavevector
                k += delta_k;
                count += 1;

                // create new decompositions for forwards and backwards analysis
                int max_wavefunction_p, max_wavefunction_m;
                DoubleHermitianEigDecomp eig_decomp_p = Diagonalise_Hamiltonian(Create_H_Using_SOIP(dft_pot, k), out max_wavefunction_p);
                DoubleHermitianEigDecomp eig_decomp_m = Diagonalise_Hamiltonian(Create_H_Using_SOIP(dft_pot, -1.0 * k), out max_wavefunction_m);

                SpinResolved_Data new_charge = new SpinResolved_Data(nx, ny);

                // cycle over each of the positive bands
                for (int i = 0; i < max_wavefunction_p; i++)
                    new_charge += Calculate_Density(eig_decomp_old_p, eig_decomp_p, i);
                // and negative bands
                for (int i = 0; i < max_wavefunction_m; i++)
                    new_charge += Calculate_Density(eig_decomp_old_m, eig_decomp_m, i);

                // check whether the density has converged
                if (Math.Abs(new_charge.Spin_Summed_Data.Min()) < eps && Math.Abs(new_charge.Spin_Summed_Data.Max()) < eps)
                    integrated = true;

                // set new eigenvalue decompositions and max_wavefunctions to the old ones
                eig_decomp_old_m = eig_decomp_m; eig_decomp_old_p = eig_decomp_p;
                max_wavefunction_old_m = max_wavefunction_m; max_wavefunction_old_p = max_wavefunction_p;

                // and finally, add the charge density calculated
                charge_density += new_charge;
            }

            kmax = k;
        }

        DoubleHermitianMatrix Create_H_Using_SOIP(Band_Data dft_pot, double k)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(2 * nx * ny);

            // create the SOIP parts of the matrix
            result = Create_SOIP_Matrix(Direction.x, Direction.y, Direction.z, k, nx, ny);
            result += Create_SOIP_Matrix(Direction.y, Direction.x, Direction.z, k, nx, ny);

            // and add the scalar terms from the k-direction
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // add momentum and potential terms
                    result[i * ny + j, i * ny + j] += Physics_Base.hbar * Physics_Base.hbar * k * k / (2.0 * Physics_Base.mass) + dft_pot.mat[i, j];
                    result[i * ny + j + nx * ny, i * ny + j + nx * ny] += Physics_Base.hbar * Physics_Base.hbar * k * k / (2.0 * Physics_Base.mass) + dft_pot.mat[i, j];

                    // add spin-orbit terms
                    result[i * ny + j, i * ny + j + nx * ny] += (-1.0 * I * alpha * dV_x.mat[i, j] - alpha * dV_y.mat[i, j]) * k;
                    result[i * ny + j + nx * ny, i * ny + j] += (I * alpha * dV_x.mat[i, j] - alpha * dV_y.mat[i, j]) * k;
                }

            return result;
        }

        private DoubleHermitianMatrix Create_SOIP_Matrix(Direction p, Direction dV, Direction sigma, double k, int nx, int ny)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(2 * nx * ny);

            double theta;
            Band_Data dV_data;

            if (dV == Direction.x)
            {
                dV_data = dV_x;
                theta = theta_x;
            }
            else if (dV == Direction.y)
            {
                dV_data = dV_y;
                theta = theta_y;
            }
            else if (dV == Direction.z)
            {
                dV_data = new Band_Data(nx, ny, 0.0);
                theta = 0.0;
            }
            else
                throw new Exception("Error - It's completely impossible to get here");


            if (p == Direction.x)
                if (sigma == Direction.x)
                    throw new InvalidArgumentException("Error - no spin-orbit interaction between p_x and sigma_x");
                else if (sigma == Direction.y)
                    throw new NotImplementedException();
                else if (sigma == Direction.z)
                {
                    for (int i = 0; i < nx; i++)
                        for(int j = 0; j < ny; j++)
                        {
                            // off-diagonal couplings in x backward
                            if (i != 0)
                            {
                                // for up -> up
                                result[i * ny + j, i * ny + j - ny] = NMathFunctions.Exp(-1.0 * I * theta * dV_data.mat[i, j]) * tx;
                                // for down -> down
                                result[i * ny + j + nx * ny, i * ny + j - ny + nx * ny] = NMathFunctions.Exp(I * theta * dV_data.mat[i, j]) * tx;
                            }
                            // off-diagonal couplings in x forward
                            if (i != nx - 1)
                            {
                                // for up -> up
                                result[i * ny + j, i * ny + j + ny] = NMathFunctions.Exp(-1.0 * I * theta * dV_data.mat[i, j]) * tx;
                                // for down -> down
                                result[i * ny + j + nx * ny, i * ny + j + ny + nx * ny] = NMathFunctions.Exp(I * theta * dV_data.mat[i, j]) * tx;
                            }
                            // and diagonal couplings
                            result[i * ny + j, i * ny + j] = -2.0 * Math.Cos(theta * dV_data.mat[i, j]) * tx;
                            result[i * ny + j + nx * ny, i * ny + j + nx * ny] = -2.0 * Math.Cos(theta * dV_data.mat[i, j]) * tx;
                        }
                }
                else
                    throw new Exception("Error - It's completely impossible to get here");
            else if (p == Direction.y)
                if (sigma == Direction.x)
                    throw new NotImplementedException();
                else if (sigma == Direction.y)
                    throw new InvalidArgumentException("Error - no spin-orbit interaction between p_y and sigma_y");
                else if (sigma == Direction.z)
                {
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                        {
                            // off-diagonal couplings in y backward
                            if (j != 0)
                            {
                                // for up -> up
                                result[i * ny + j, i * ny + j - 1] = NMathFunctions.Exp(-1.0 * I * theta * dV_data.mat[i, j]) * ty;
                                // for down -> down
                                result[i * ny + j + nx * ny, i * ny + j - 1 + nx * ny] = NMathFunctions.Exp(I * theta * dV_data.mat[i, j]) * ty;
                            }
                            // off-diagonal couplings in y forward
                            if (j != ny - 1)
                            {
                                // for up -> up
                                result[i * ny + j, i * ny + j + 1] = NMathFunctions.Exp(-1.0 * I * theta * dV_data.mat[i, j]) * ty;
                                // for down -> down
                                result[i * ny + j + nx * ny, i * ny + j + 1 + nx * ny] = NMathFunctions.Exp(I * theta * dV_data.mat[i, j]) * ty;
                            }
                            // and diagonal couplings
                            result[i * ny + j, i * ny + j] = -2.0 * Math.Cos(theta * dV_data.mat[i, j]) * ty;
                            result[i * ny + j + nx * ny, i * ny + j + nx * ny] = -2.0 * Math.Cos(theta * dV_data.mat[i, j]) * ty;
                        }
                }
                else
                    throw new Exception("Error - It's completely impossible to get here");
            else if (p == Direction.z)
                throw new NotImplementedException("Error - at the moment, when momentum is in the direction of the wire, you should use a different method");
            else
                throw new Exception("Error - It's completely impossible to get here");

            return result;
        }

        DoubleHermitianMatrix Create_Hamiltonian(Band_Data dft_pot, double k)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(2 * nx * ny);
            DoubleHermitianMatrix h_upup, h_updn, h_dnup, h_dndn;

            // create sub-matrices with no spin orbit
            h_upup = Create_NoSOI_Hamiltonian(dft_pot, k, nx, ny);      // takes up to up (ie. H_11)
            h_dndn = Create_NoSOI_Hamiltonian(dft_pot, k, nx, ny);      // takes dn to dn (ie. H_22)
            h_updn = new DoubleHermitianMatrix(nx * ny);                                // takes up to dn (ie. H_21)
            h_dnup = new DoubleHermitianMatrix(nx * ny);                                // takes dn to up (ie. H_12)

            // add spin-orbit
            h_upup -= alpha * Create_SOI_Hamiltonian(Direction.y, Direction.x, k, nx, ny);
            h_dndn -= -1.0 * alpha * Create_SOI_Hamiltonian(Direction.y, Direction.x, k, nx, ny);
            h_updn -= alpha * (Create_SOI_Hamiltonian(Direction.z, Direction.y, nx, ny) - Create_SOI_Hamiltonian(Direction.y, Direction.z, nx, ny) 
                                + I * Create_SOI_Hamiltonian(Direction.z, Direction.x, k, nx, ny));
            h_dnup -= alpha * (Create_SOI_Hamiltonian(Direction.z, Direction.y, nx, ny) - Create_SOI_Hamiltonian(Direction.y, Direction.z, nx, ny)
                                - I * Create_SOI_Hamiltonian(Direction.z, Direction.x, k, nx, ny));

            // recombine matrices and output result
            for (int i = 0; i < nx * ny; i++)
                for (int j = 0; j < nx * ny; j++)
                {
                    result[i, j] = h_upup[i, j];
                    result[i, j + nx * ny] = h_dnup[i, j];
                    result[i + nx * ny, j] = h_updn[i, j];
                    result[i + nx * ny, j + nx * ny] = h_dndn[i, j];
                }

            return result;
        }

        DoubleHermitianMatrix Create_SOI_Hamiltonian(Direction dV, Direction p, double k, int nx, int ny)
        {
            DoubleHermitianMatrix result = 0.0 * new DoubleHermitianMatrix(nx * ny);
            Band_Data dV_data;

            if (dV == Direction.x)
                return result;
            else if (dV == Direction.y)
                dV_data = dV_x;
            else if (dV == Direction.z)
                dV_data = dV_y;
            else
                throw new Exception("Error - It's completely impossible to get here");

            if (p == Direction.x)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result[i * ny + j, i * ny + j] = Physics_Base.hbar * k * dV_data.mat[i, j];
            else if (p == Direction.y)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        if (i != nx - 1)
                            result[i * ny + j, i * ny + j + ny] = -0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dx;
                        if (i != 0)
                            result[i * ny + j, i * ny + j - ny] = 0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dx;
                    }
            else if (p == Direction.z)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        if (i != nx - 1)
                            result[i * ny + j, i * ny + j + 1] = -0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dy;
                        if (i != 0)
                            result[i * ny + j, i * ny + j - 1] = 0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dy;
                    }
            else
                throw new Exception("Error - It's completely impossible to get here");

            return result;
        }

        DoubleHermitianMatrix Create_SOI_Hamiltonian(Direction dV, Direction p, int nx, int ny)
        {
            return Create_SOI_Hamiltonian(dV, p, 0.0, nx, ny);
        }

        DoubleHermitianMatrix Create_NoSOI_Hamiltonian(Band_Data dft_pot, double k, int nx, int ny)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(nx * ny);

            // set off diagonal elements 
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // coupling sites in the growth direction
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

            double[,] potential = new double[nx, ny];
            // set diagonal elements
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    potential[i, j] = dft_pot.mat[i, j];
                    potential[i, j] += Physics_Base.hbar * Physics_Base.hbar * k * k / (2.0 * Physics_Base.mass);
                    result[i * ny + j, i * ny + j] = -2.0 * tx + -2.0 * ty + potential[i, j];
                }

            return result;
        }

        bool Check_Integration_Convergence(SpinResolved_Data density_k, double eps)
        {
            // checks whether the maximum (negative) density is less than eps
            double max_dens = 0.0;
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                    if (density_k.Spin_Up.mat[i, j] + density_k.Spin_Down.mat[i, j] < max_dens)
                        max_dens = density_k.Spin_Up.mat[i, j] + density_k.Spin_Down.mat[i, j];

            return eps > -1.0 * max_dens;
        }

        SpinResolved_Data Calculate_Density(DoubleHermitianEigDecomp eig_decomp_old, DoubleHermitianEigDecomp eig_decomp_new, int wavefunction)
        {
                DoubleMatrix dens_up = new DoubleMatrix(nx, ny, 0.0);
                DoubleMatrix dens_down = new DoubleMatrix(nx, ny, 0.0);

                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        // do not add anything to the density if on the edge of the domain
                        if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                            continue;

                        // calculate the density of this spin configuration this hamiltonian
                        dens_up[i, j] += g_1D * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j]) * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j]) * Interpolate_Fermi_Function(eig_decomp_old, eig_decomp_new, wavefunction);
                        dens_down[i, j] += g_1D * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j + nx * ny]) * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j + nx * ny]) * Interpolate_Fermi_Function(eig_decomp_old, eig_decomp_new, wavefunction);
                    }

            // and multiply the density by -e to get the charge density (as these are electrons)
            return -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));;
        }

        double Interpolate_Fermi_Function(DoubleHermitianEigDecomp eig_decomp_old, DoubleHermitianEigDecomp eig_decomp_new, int wavefunction)
        {
            double E0 = eig_decomp_old.EigenValue(wavefunction);
            double dE_dk = (eig_decomp_new.EigenValue(wavefunction) - eig_decomp_old.EigenValue(wavefunction)) / delta_k;

            double beta = 1.0 / (Physics_Base.kB * temperature);
            OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double k) => Math.Pow(Math.Exp(beta * (E0 + dE_dk * k)) + 1, -1.0)));
            dos_integrand.Integrator = new GaussKronrodIntegrator();

            return dos_integrand.Integrate(0, delta_k);
        }

        DoubleHermitianEigDecomp Diagonalise_Hamiltonian(DoubleHermitianMatrix hamiltonian, out int max_wavefunction)
        {
            DoubleHermitianEigDecomp eig_decomp;

            DoubleHermitianEigDecompServer eig_server = new DoubleHermitianEigDecompServer();
            eig_server.ComputeEigenValueRange(E_min, no_kb_T * Physics_Base.kB * temperature);
            eig_server.ComputeVectors = true;
            eig_decomp = eig_server.Factor(hamiltonian);

            max_wavefunction = 0;
            if (eig_decomp.EigenValues.Length != 0)
            {
                double min_eigval = eig_decomp.EigenValues.Min();
                max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;
            }

            return eig_decomp;
        }

        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data charge_density_deriv, Band_Data chem_pot)
        {
            Get_SOI_parameters(chem_pot);

            //   if (dV_x == null)
            //       throw new Exception("Error - Band structure derivatives are null!  Have you initiated this type properly by calling Get_SOI_parameters(Band_Data chem_pot)?");

            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            // reset charge density
            charge_density_deriv = 0.0 * charge_density_deriv;

            // integrate up density from k = 0 to k = k_F
            int count = 0;
            double k = 0;

            // create initial hamiltonian and eigendecomposition
            int max_wavefunction_old_p, max_wavefunction_old_m;
            DoubleHermitianMatrix hamiltonian_p = Create_H_Using_SOIP(dft_pot, k);
            DoubleHermitianEigDecomp eig_decomp_old_p = Diagonalise_Hamiltonian(hamiltonian_p, out max_wavefunction_old_p);
            DoubleHermitianEigDecomp eig_decomp_old_m = Diagonalise_Hamiltonian(hamiltonian_p, out max_wavefunction_old_m);

            while (k < kmax)
            {
                // increment the wavevector
                k += delta_k;
                count += 1;

                // create new decompositions for forwards and backwards analysis
                int max_wavefunction_p, max_wavefunction_m;
                DoubleHermitianEigDecomp eig_decomp_p = Diagonalise_Hamiltonian(Create_H_Using_SOIP(dft_pot, k), out max_wavefunction_p);
                DoubleHermitianEigDecomp eig_decomp_m = Diagonalise_Hamiltonian(Create_H_Using_SOIP(dft_pot, -1.0 * k), out max_wavefunction_m);

                SpinResolved_Data new_charge_deriv = new SpinResolved_Data(nx, ny);

                // cycle over each of the positive bands
                for (int i = 0; i < max_wavefunction_p; i++)
                    new_charge_deriv += Calculate_Density_Derivative(eig_decomp_old_p, eig_decomp_p, i);
                // and negative bands
                for (int i = 0; i < max_wavefunction_m; i++)
                    new_charge_deriv += Calculate_Density_Derivative(eig_decomp_old_m, eig_decomp_m, i);

                // set new eigenvalue decompositions and max_wavefunctions to the old ones
                eig_decomp_old_m = eig_decomp_m; eig_decomp_old_p = eig_decomp_p;
                max_wavefunction_old_m = max_wavefunction_m; max_wavefunction_old_p = max_wavefunction_p;

                // and finally, add the charge density calculated
                charge_density_deriv += new_charge_deriv;
            }
        }

        SpinResolved_Data Calculate_Density_Derivative(DoubleHermitianEigDecomp eig_decomp_old, DoubleHermitianEigDecomp eig_decomp_new, int wavefunction)
        {
            DoubleMatrix dens_up_deriv = new DoubleMatrix(nx, ny, 0.0);
            DoubleMatrix dens_down_deriv = new DoubleMatrix(nx, ny, 0.0);

            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    // calculate the density of this spin configuration this hamiltonian
                    dens_up_deriv[i, j] += g_1D * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j]) * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j]) * Interpolate_Fermi_Function_Derivative(eig_decomp_old, eig_decomp_new, wavefunction);
                    dens_down_deriv[i, j] += g_1D * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j + nx * ny]) * DoubleComplex.Norm(eig_decomp_old.EigenVector(wavefunction)[i * ny + j + nx * ny]) * Interpolate_Fermi_Function_Derivative(eig_decomp_old, eig_decomp_new, wavefunction);
                }

            // and multiply the density by -e to get the charge density (as these are electrons)
            return Physics_Base.q_e * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up_deriv), new Band_Data(dens_down_deriv));
        }

        double Interpolate_Fermi_Function_Derivative(DoubleHermitianEigDecomp eig_decomp_old, DoubleHermitianEigDecomp eig_decomp_new, int wavefunction)
        {
            double E0 = eig_decomp_old.EigenValue(wavefunction);
            double dE_dk = (eig_decomp_new.EigenValue(wavefunction) - eig_decomp_old.EigenValue(wavefunction)) / delta_k;

            double beta = 1.0 / (Physics_Base.kB * temperature);
            OneVariableFunction dos_integrand = new OneVariableFunction((Func<double, double>)((double k) => -1.0 * beta * Math.Exp(beta * (E0 + dE_dk * k)) * Math.Pow(Math.Exp(beta * (E0 + dE_dk * k)) + 1, -2.0)));
            dos_integrand.Integrator = new GaussKronrodIntegrator();

            return dos_integrand.Integrate(0, delta_k);
        }

        /// <summary>
        /// Calculates and prints the band structure
        /// </summary>
        void Print_Band_Structure(Band_Data dft_pot, ILayer[] layers, int Nk, double dk, string outfile, int max_eigval)
        {
            // set temperature high so that all wave functions will be calculated
            double old_temperature = temperature;
            temperature = 1000.0;

            if (dV_x == null)
                throw new Exception("Error - Band structure derivatives are null!  Have you initiated this type properly by calling Get_SOI_parameters(Band_Data chem_pot)?");

            // calculate the energies up to a given maximum energy of the lowest state
            double k = 0;
            double[][] energies = new double[Nk][];
            for (int i = 0; i < Nk; i++)
            {
                k = i * dk;
                Console.WriteLine(i.ToString() + ": Calculating for k = " + k.ToString());
                // generate the Hamiltonian for this k value
                DoubleHermitianMatrix hamiltonian = Create_H_Using_SOIP(dft_pot, k);
                // and diagonalise it
                int max_wavefunction;
                DoubleHermitianEigDecomp eig_decomp = Diagonalise_Hamiltonian(hamiltonian, out max_wavefunction);

                // add the calculated energies up to either the maximum required eigenvalue or
                // to the maximum calculated wave function (which is 50*kb*T above the chemical potential)
                double[] tmp_energies = new double[max_eigval];
                for (int j = 0; j < max_eigval; j++)
                    if (j < max_wavefunction)
                        tmp_energies[j] = eig_decomp.EigenValues[j];// - 0.5 * Physics_Base.hbar * Physics_Base.hbar * k * k / Physics_Base.mass;
                    else
                        tmp_energies[j] = eig_decomp.EigenValues[max_wavefunction - 1];// -0.5 * Physics_Base.hbar * Physics_Base.hbar * k * k / Physics_Base.mass;

                energies[i] = tmp_energies;
            }

            // output the data to file
            StreamWriter sw = new StreamWriter(outfile);
            for (int i = 0; i < Nk; i++)
            {
                for (int j = 0; j < max_eigval; j++)
                    sw.Write(energies[i][j].ToString() + "\t");

                sw.WriteLine();
            }
            sw.Close();

            // reset temperature
            temperature = old_temperature;
        }

        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot)
        {
            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            // diagonalise matrix for k = 0
            int tmp = 0;
            DoubleHermitianEigDecomp eig_decomp = Diagonalise_Hamiltonian(Create_Hamiltonian(dft_pot, 0.0), out tmp);

            return eig_decomp.EigenValues;
        }
    }
}
