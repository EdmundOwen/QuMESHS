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
        double delta_k = 0.03;           // integration variable in k for band structure
        double eps = 1.0e-9;            // minimum density added before convergence

        DoubleComplex I = new DoubleComplex(0.0, 1.0);

        double ty, tz;
        double alpha;
        double h;
        double g_1D;
        double E_min = -10.0;           // minimum energy value for Lanczos matrix diagonalisation routine in meV
                                        // should be less than minimum energy expected for any k-vector

        Band_Data dV_y, dV_z, dV_yz;

        public TwoD_SO_DFTSolver(Experiment exp)
            : base(exp)
        {
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);
            double r_so = 117.1;
            alpha = r_so / (Physics_Base.q_e * Physics_Base.hbar);

            h = 3.0 * delta_k / 8.0;
            g_1D = 2.0 / Math.PI;
        }

        public void Get_SOI_parameters(Band_Data chem_pot)
        {
            // initialise matrices with same dimension as the chem_pot data which we will receive later
            dV_y = new Band_Data(new DoubleMatrix(chem_pot.mat.Rows, chem_pot.mat.Cols));
            dV_z = new Band_Data(new DoubleMatrix(chem_pot.mat.Rows, chem_pot.mat.Cols));
            dV_yz = new Band_Data(new DoubleMatrix(chem_pot.mat.Rows, chem_pot.mat.Cols));

            // calculate the derivatives of chem_pot in y and z and the second mixed derivative
            for (int i = 1; i < ny - 1; i++)
                for (int j = 1; j < nz - 1; j++)
                {
                    dV_y.mat[i, j] = (chem_pot.mat[i + 1, j] - chem_pot.mat[i - 1, j]) / (2.0 * dy);
    //                dV_z.mat[i, j] = (chem_pot.mat[i, j + 1] - chem_pot.mat[i, j - 1]) / (2.0 * dz);
                    dV_yz.mat[i, j] = (chem_pot.mat[i + 1, j + 1] - chem_pot.mat[i - 1, j + 1] - chem_pot.mat[i + 1, j - 1] + chem_pot.mat[i - 1, j - 1]) / (4.0 * dy * dz);
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

            // test
            Console.WriteLine("TEST!!!!!!!!!!!!!!!!!!!!!");
            
            dy = 0.1;
            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            dz = 0.1;
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);

            double l_w = ny * dy / 10.0;
            double omega_z = 50000.0;
            double omega = Physics_Base.hbar / Physics_Base.mass / l_w / l_w;
            for (int i = 0; i < dft_pot.mat.Rows; i++)
                for (int j = 0; j < dft_pot.mat.Cols; j++)
                    dft_pot.mat[i, j] = 0.5 * Physics_Base.mass * omega * omega * (i - (ny + 1) / 2) * (i - (ny + 1) / 2) * dy * dy + 0.5 * Physics_Base.mass * omega_z * omega_z * (j - (nz + 1) / 2) * (j - (nz + 1) / 2) * dz * dz;

            alpha = 1.0;
            double beta = Physics_Base.hbar * Physics_Base.hbar / (2.0 * Physics_Base.mass * 10.0 * l_w);
            double factor = beta / (Physics_Base.hbar * l_w * Physics_Base.mass * omega * omega);
            Get_SOI_parameters(factor * dft_pot);

            temperature = 1000.0;

            Console.WriteLine("hbar omega / 2 = " + (Physics_Base.hbar * omega * 0.5).ToString());
            Print_Band_Structure(dft_pot, layers, 21, 1.0 / l_w, "tmp.dat", 20);
            Console.WriteLine("End of test");

            // reset charge density
            charge_density = 0.0 * charge_density;
            
            // integrate up density from k = 0 to k = k_F
            bool integrated = false;
            int count = 0;
            double k = 0;
            while (!integrated)
            {
                // generate the Hamiltonian for this k value
                DoubleHermitianMatrix hamiltonian_p = Create_Hamiltonian(layers, dft_pot, k);
                DoubleHermitianMatrix hamiltonian_m = Create_Hamiltonian(layers, dft_pot, -1.0 * k);

                // calculate density for this value of k
                SpinResolved_Data charge_density_k = Calculate_Density(hamiltonian_p);
                if (k != 0)
                    charge_density_k += Calculate_Density(hamiltonian_m);

                // add new k density according to Simpson's 3/8th rule
                if (count == 0)
                {
                    charge_density = h * charge_density_k;
                    // and check to see whether there is any density being added
                    integrated = Check_Convergence(charge_density_k, eps);
                }
                else if (count % 3 != 0)
                    charge_density += 3.0 * h * charge_density_k;
                else
                {
                    // check whether there was any density at this wave vector to add
                    integrated = Check_Integration_Convergence(charge_density_k, eps);

                    if (integrated)
                        charge_density += h * charge_density_k;
                    else
                        charge_density += 2.0 * h * charge_density_k;
                }

                k += delta_k;
                count += 1;
            }
        }

        DoubleHermitianMatrix Create_Hamiltonian(ILayer[] layers, Band_Data dft_pot, double k)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(2 * ny * nz);
            DoubleHermitianMatrix h_upup, h_updn, h_dnup, h_dndn;

            // create sub-matrices with no spin orbit
            h_upup = Create_NoSOI_Hamiltonian(dft_pot, k, ny, nz);      // takes up to up (ie. H_11)
            h_dndn = Create_NoSOI_Hamiltonian(dft_pot, k, ny, nz);      // takes dn to dn (ie. H_22)
            h_updn = new DoubleHermitianMatrix(ny * nz);                                // takes up to dn (ie. H_21)
            h_dnup = new DoubleHermitianMatrix(ny * nz);                                // takes dn to up (ie. H_12)

            // add spin-orbit
            h_upup -= alpha * Create_SOI_Hamiltonian(Direction.y, Direction.x, k, ny, nz);
            h_dndn -= -1.0 * alpha * Create_SOI_Hamiltonian(Direction.y, Direction.x, k, ny, nz);
            h_updn -= alpha * (Create_SOI_Hamiltonian(Direction.z, Direction.y, ny, nz) - Create_SOI_Hamiltonian(Direction.y, Direction.z, ny, nz) 
                                + I * Create_SOI_Hamiltonian(Direction.z, Direction.x, k, ny, nz));
            h_dnup -= alpha * (Create_SOI_Hamiltonian(Direction.z, Direction.y, ny, nz) - Create_SOI_Hamiltonian(Direction.y, Direction.z, ny, nz)
                                - I * Create_SOI_Hamiltonian(Direction.z, Direction.x, k, ny, nz));

            // recombine matrices and output result
            for (int i = 0; i < ny * nz; i++)
                for (int j = 0; j < ny * nz; j++)
                {
                    result[i, j] = h_upup[i, j];
                    result[i, j + ny * nz] = h_dnup[i, j];
                    result[i + ny * nz, j] = h_updn[i, j];
                    result[i + ny * nz, j + ny * nz] = h_dndn[i, j];
                }

            return result;
        }

        DoubleHermitianMatrix Create_SOI_Hamiltonian(Direction dV, Direction p, double k, int ny, int nz)
        {
            DoubleHermitianMatrix result = 0.0 * new DoubleHermitianMatrix(ny * nz);
            Band_Data dV_data;

            if (dV == Direction.x)
                return result;
            else if (dV == Direction.y)
                dV_data = dV_y;
            else if (dV == Direction.z)
                dV_data = dV_z;
            else
                throw new Exception("Error - It's completely impossible to get here");

            if (p == Direction.x)
                for (int i = 0; i < ny; i++)
                    for (int j = 0; j < nz; j++)
                        result[i * nz + j, i * nz + j] = Physics_Base.hbar * k * dV_data.mat[i, j];
            else if (p == Direction.y)
                for (int i = 0; i < ny; i++)
                    for (int j = 0; j < nz; j++)
                    {
                        if (i != ny - 1)
                            result[i * nz + j, i * nz + j + nz] = -0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dy;
                        if (i != 0)
                            result[i * nz + j, i * nz + j - nz] = 0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dy;
                    }
            else if (p == Direction.z)
                for (int i = 0; i < ny; i++)
                    for (int j = 0; j < nz; j++)
                    {
                        if (i != ny - 1)
                            result[i * nz + j, i * nz + j + 1] = -0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dz;
                        if (i != 0)
                            result[i * nz + j, i * nz + j - 1] = 0.5 * I * Physics_Base.hbar * dV_data.mat[i, j] / dz;
                    }
            else
                throw new Exception("Error - It's completely impossible to get here");

            return result;
        }

        DoubleHermitianMatrix Create_SOI_Hamiltonian(Direction dV, Direction p, int ny, int nz)
        {
            return Create_SOI_Hamiltonian(dV, p, 0.0, ny, nz);
        }

        DoubleHermitianMatrix Create_NoSOI_Hamiltonian(Band_Data dft_pot, double k, int ny, int nz)
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
                    potential[i, j] = dft_pot.mat[i, j];
                    potential[i, j] += Physics_Base.hbar * Physics_Base.hbar * k * k / (2.0 * Physics_Base.mass);
                    result[i * nz + j, i * nz + j] = -2.0 * ty + -2.0 * tz + potential[i, j];
                }

            return result;
        }

        bool Check_Integration_Convergence(SpinResolved_Data density_k, double eps)
        {
            // checks whether the maximum (negative) density is less than eps
            double max_dens = 0.0;
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    if (density_k.Spin_Up.mat[i, j] + density_k.Spin_Down.mat[i, j] < max_dens)
                        max_dens = density_k.Spin_Up.mat[i, j] + density_k.Spin_Down.mat[i, j];

            return eps > -1.0 * max_dens;
        }

        SpinResolved_Data Calculate_Density(DoubleHermitianMatrix hamiltonian)
        {
            int max_wavefunction;
            DoubleHermitianEigDecomp eig_decomp = Diagonalise_Hamiltonian(hamiltonian, out max_wavefunction);

            DoubleMatrix dens_up = new DoubleMatrix(ny, nz, 0.0);
            DoubleMatrix dens_down = new DoubleMatrix(ny, nz, 0.0);

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                        continue;

                    for (int n = 0; n < max_wavefunction; n++)
                    {
                        // calculate the density of this spin configuration this hamiltonian
                        dens_up[i, j] += g_1D * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * nz + j]) * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * nz + j]) * Physics_Base.Get_Fermi_Function(eig_decomp.EigenValue(n), 0.0, temperature);
                        dens_down[i, j] += g_1D * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * nz + j + ny * nz]) * DoubleComplex.Norm(eig_decomp.EigenVector(n)[i * nz + j + ny * nz]) * Physics_Base.Get_Fermi_Function(eig_decomp.EigenValue(n), 0.0, temperature);
                    }
                }

            // and multiply the density by -e to get the charge density (as these are electrons)
            return -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        DoubleHermitianEigDecomp Diagonalise_Hamiltonian(DoubleHermitianMatrix hamiltonian, out int max_wavefunction)
        {
            DoubleHermitianEigDecomp eig_decomp;

            DoubleHermitianEigDecompServer eig_server = new DoubleHermitianEigDecompServer();
  //          eig_server.ComputeEigenValueRange(E_min, no_kb_T * Physics_Base.kB * temperature);
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

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    // leave the edges zeroed
                    if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                        continue;

                    double y = dy * i + ymin;
                    double z = dz * j + zmin;

                    // get the relevant layer and if it's frozen out, don't recalculate the dopent charge
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, y, z);

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

        /// <summary>
        /// Calculates and prints the band structure
        /// </summary>
        void Print_Band_Structure(Band_Data dft_pot, ILayer[] layers, int Nk, double dk, string outfile, int max_eigval)
        {
            if (dV_y == null)
                throw new Exception("Error - Band structure derivatives are null!  Have you initiated this type properly by calling Get_SOI_parameters(Band_Data chem_pot)?");

            // calculate the energies up to a given maximum energy of the lowest state
            double k = 0;
            double[][] energies = new double[Nk][];
            for (int i = 0; i < Nk; i++)
            {
                k = i * dk;
                Console.WriteLine(i.ToString() + ": Calculating for k = " + k.ToString());
                // generate the Hamiltonian for this k value
                DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, dft_pot, k);
                // and diagonalise it
                int max_wavefunction;
                DoubleHermitianEigDecomp eig_decomp = Diagonalise_Hamiltonian(hamiltonian, out max_wavefunction);
                max_wavefunction = max_eigval;

                // add the calculated energies up to either the maximum required eigenvalue or
                // to the maximum calculated wave function (which is 50*kb*T above the chemical potential)
                double[] tmp_energies = new double[max_eigval];
                for (int j = 0; j < max_eigval; j++)
                    if (j < max_wavefunction)
                        tmp_energies[j] = eig_decomp.EigenValues[j] - 0.5 * Physics_Base.hbar * Physics_Base.hbar * k * k / Physics_Base.mass;
                    else
                        tmp_energies[j] = eig_decomp.EigenValues[max_wavefunction - 1] -0.5 * Physics_Base.hbar * Physics_Base.hbar * k * k / Physics_Base.mass;

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
        }

    }
}
