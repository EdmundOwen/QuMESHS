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
    class TwoD_SO_DFTSolver : TwoD_Density_Base
    {
        double no_kb_T = 50.0;          // number of kb_T to integrate to
        double delta_k = 0.1;           // integration variable in k for band structure
        double eps = 1.0e-9;            // minimum density added before convergence

        DoubleComplex I = new DoubleComplex(0.0, 1.0);

        double tx, ty;
        double alpha;
        double h;
        double g_1D;
        double E_min = -10.0;           // minimum energy value for Lanczos matrix diagonalisation routine in meV
                                        // should be less than minimum energy expected for any k-vector

        Band_Data dV_x, dV_y, dV_xy;

        public TwoD_SO_DFTSolver(IExperiment exp)
            : base(exp)
        {
            tx = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            double r_so = 1.0e-8;
            alpha = r_so * Physics_Base.hbar * Physics_Base.hbar / (4.0 * Physics_Base.mass * Physics_Base.mass);

            h = 3.0 * delta_k / 8.0;
            g_1D = 2.0 / Math.PI;
        }

        public void Get_SOI_parameters(Band_Data chem_pot)
        {
            // initialise matrices with same dimension as the chem_pot data which we will receive later
            dV_x = new Band_Data(new DoubleMatrix(chem_pot.mat.Rows, chem_pot.mat.Cols));
            dV_y = new Band_Data(new DoubleMatrix(chem_pot.mat.Rows, chem_pot.mat.Cols));
            dV_xy = new Band_Data(new DoubleMatrix(chem_pot.mat.Rows, chem_pot.mat.Cols));

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
            if (dV_y == null)
                throw new Exception("Error - Band structure derivatives are null!  Have you initiated this type properly by calling Get_SOI_parameters(Band_Data chem_pot)?");

            // convert the chemical potential into a quantum mechanical potential
            Band_Data dft_pot = chem_pot.DeepenThisCopy();
            Get_Potential(ref dft_pot, layers);

            // reset charge density
            charge_density = 0.0 * charge_density;
            
            // integrate up density from k = 0 to k = k_F
            bool integrated = false;
            int count = 0;
            double k = 0;
            while (!integrated)
            {
                // generate the Hamiltonian for this k value
                DoubleHermitianMatrix hamiltonian_p = Create_Hamiltonian(layers, charge_density, dft_pot, k);
                DoubleHermitianMatrix hamiltonian_m = Create_Hamiltonian(layers, charge_density, dft_pot, -1.0 * k);

                // calculate density for this value of k
                SpinResolved_Data density_k_p = Calculate_Density(hamiltonian_p);
                SpinResolved_Data density_k_m = Calculate_Density(hamiltonian_m);
                SpinResolved_Data density_k = density_k_p + density_k_m;

                // add new k density according to Simpson's 3/8th rule
                if (count % 3 != 0)
                    charge_density += 3.0 * h * density_k;
                else if (count == 0)
                    charge_density = h * density_k;
                else
                {
                    // check whether there was any density at this wave vector to add
                    integrated = Check_Integration_Convergence(density_k, eps);

                    if (integrated)
                        charge_density = h * density_k;
                    else
                        charge_density = 2.0 * h * density_k;
                }

                k += delta_k;
                count += 1;
            }
        }

        DoubleHermitianMatrix Create_Hamiltonian(ILayer[] layers, SpinResolved_Data charge_density, Band_Data dft_pot, double k)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(2 * nx * ny);
            DoubleHermitianMatrix h_upup, h_updn, h_dnup, h_dndn;

            // create sub-matrices with no spin orbit
            h_upup = Create_NoSOI_Hamiltonian(dft_pot, charge_density, k, nx, ny);      // takes up to up (ie. H_11)
            h_dndn = Create_NoSOI_Hamiltonian(dft_pot, charge_density, k, nx, ny);      // takes dn to dn (ie. H_22)
            h_updn = new DoubleHermitianMatrix(nx * ny);                                // takes up to dn (ie. H_21)
            h_dnup = new DoubleHermitianMatrix(nx * ny);                                // takes dn to up (ie. H_12)

            // add spin-orbit
            h_upup += alpha * Create_SOI_Hamiltonian(Direction.x, Direction.z, k, nx, ny);
            h_dndn += -1.0 * alpha * Create_SOI_Hamiltonian(Direction.x, Direction.z, k, nx, ny);
            h_updn += alpha * (Create_SOI_Hamiltonian(Direction.y, Direction.x, nx, ny) - Create_SOI_Hamiltonian(Direction.x, Direction.y, nx, ny) 
                                + I * Create_SOI_Hamiltonian(Direction.y, Direction.z, k, nx, ny));
            h_dnup += alpha * (Create_SOI_Hamiltonian(Direction.y, Direction.x, nx, ny) - Create_SOI_Hamiltonian(Direction.x, Direction.y, nx, ny)
                                - I * Create_SOI_Hamiltonian(Direction.y, Direction.z, k, nx, ny));

            // recombine matrices and output result
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
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

            if (dV == Direction.z)
                return result;
            else if (dV == Direction.x)
                dV_data = dV_x;
            else if (dV == Direction.y)
                dV_data = dV_y;
            else
                throw new Exception("Error - It's completely impossible to get here");

            if (p == Direction.z)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result[i, j] = k * dV_data.mat[i, j];
            else if (p == Direction.x)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        if (i != nx - 1)
                            result[i * ny + j, i * ny + j + ny] = 0.5 * dV_data.mat[i, j] / dx;
                        if (i != 0)
                            result[i * ny + j, i * ny + j - ny] = -0.5 * dV_data.mat[i, j] / dx;
                    }
            else if (p == Direction.y)
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                    {
                        if (i != nx - 1)
                            result[i * ny + j, i * ny + j + 1] = 0.5 * dV_data.mat[i, j] / dx;
                        if (i != 0)
                            result[i * ny + j, i * ny + j - 1] = -0.5 * dV_data.mat[i, j] / dx;
                    }
            else
                throw new Exception("Error - It's completely impossible to get here");

            return result;
        }

        DoubleHermitianMatrix Create_SOI_Hamiltonian(Direction dV, Direction p, int nx, int ny)
        {
            return Create_SOI_Hamiltonian(dV, p, 0.0, nx, ny);
        }

        DoubleHermitianMatrix Create_NoSOI_Hamiltonian(Band_Data dft_pot, SpinResolved_Data charge_density, double k, int nx, int ny)
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
                    potential[i, j] = dft_pot.mat[i, j];// +Physics_Base.Get_XC_Potential(charge_density.Spin_Summed_Data.mat[i, j]);  This should already be included in the input chemical potential
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

        SpinResolved_Data Calculate_Density(DoubleHermitianMatrix hamiltonian)
        {
            DoubleHermitianEigDecompServer eig_server = new DoubleHermitianEigDecompServer();
            eig_server.ComputeEigenValueRange(E_min, no_kb_T * Physics_Base.kB * temperature);
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
                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    for (int n = 0; n < max_wavefunction; n++)
                    {
                        // calculate the density of this spin configuration this hamiltonian
                        dens_up[i, j] += g_1D * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * ny + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * ny + j]) * Physics_Base.Get_Fermi_Function(eig_decomp.EigenValue(n), 0.0, temperature);
                        dens_down[i, j] += g_1D * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * ny + j + nx * ny]) * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * ny + j + nx * ny]) * Physics_Base.Get_Fermi_Function(eig_decomp.EigenValue(n), 0.0, temperature);
                    }
                }

            // and multiply the density by -e to get the charge density (as these are electrons)
            return -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // leave the edges zeroed
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    double x = dx * i + xmin;
                    double y = dy * j + ymin;

                    // get the relevant layer and if it's frozen out, don't recalculate the dopent charge
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, plane, x, y, pos_z);

                    ZeroD_Density charge_calc = new ZeroD_Density(current_Layer, temperature);
                    if (!current_Layer.Dopents_Frozen_Out(temperature))
                    {
                        double local_dopent_density_deriv = charge_calc.Get_DopentDensityDeriv(chem_pot.mat[i, j]);
                        dopent_density_deriv.Spin_Up.mat[i, j] = 0.5 * local_dopent_density_deriv;
                        dopent_density_deriv.Spin_Down.mat[i, j] = 0.5 * local_dopent_density_deriv;
                    }
                    else
                    {
                        dopent_density_deriv.Spin_Up.mat[i, j] = 0.0;
                        dopent_density_deriv.Spin_Down.mat[i, j] = 0.0;
                    }

                    carrier_density_deriv.Spin_Up.mat[i, j] = 0.5 * charge_calc.Get_CarrierDensityDeriv(chem_pot.mat[i, j]);
                    carrier_density_deriv.Spin_Down.mat[i, j] = 0.5 * charge_calc.Get_CarrierDensityDeriv(chem_pot.mat[i, j]);
                }

            return carrier_density_deriv + dopent_density_deriv;
        }
        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }
    }
}
