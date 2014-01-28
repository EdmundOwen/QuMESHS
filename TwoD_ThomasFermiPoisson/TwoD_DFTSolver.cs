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
    class TwoD_DFTSolver : Density_Base
    {
        Experiment exp;

        double dos_cutoff = 100;    // density of states cutoff to prevent divergences at the band edges
        double no_kb_T = 20;        // number of kb_T to integrate to

        int tmp_yval, tmp_zval;
        Double tmp_eigval;
        DoubleComplexVector tmp_eigvec;

        double ty, tz;

        public TwoD_DFTSolver(Experiment exp)
            : base(exp.Temperature, exp.Dx_Dens, exp.Dy_Dens, exp.Dz_Dens, exp.Nx_Dens, exp.Ny_Dens, exp.Nz_Dens, exp.Xmin_Dens, exp.Ymin_Dens, exp.Zmin_Dens)
        {
            this.exp = exp;

            ty = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
            tz = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data pot)
        {
            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, density, pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            int max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            DoubleMatrix dens_up = new DoubleMatrix(ny, nz, 0.0);
            DoubleMatrix dens_down = new DoubleMatrix(ny, nz, 0.0);
            OneVariableFunction dens_of_states = new OneVariableFunction(new Func<double, double>(Density_Of_States));

            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    double dens_val = 0.0;

                    for (int k = 0; k < max_wavefunction; k++)
                    {
                        // set temporary eigenvalue and eigenvector
                        tmp_eigval = eig_decomp.EigenValue(k); tmp_eigvec = eig_decomp.EigenVector(k);
                        // and position
                        tmp_yval = i;
                        tmp_zval = j;

                        // and integrate the density of states at this position for this eigenvector from the minimum energy to
                        // (by default) 20 * k_b * T above mu = 0
                        dens_val += dens_of_states.Integrate(min_eigval, no_kb_T * Physics_Base.kB * temperature);
                    }

                    // just share the densities (there is no spin polarisation)
                    dens_up[i, j] = 0.5 * dens_val;
                    dens_down[i, j] = 0.5 * dens_val;
                }

            //
            //Console.Clear();
            //for (int i = 0; i < nz; i++)
            //    Console.WriteLine(charge_density.Spin_Summed_Data[i].ToString() + "\t" + (-1.0 * Physics_Base.q_e * (dens_up[i] + dens_down[i])).ToString() + "\t" + (chem_pot.vec[i] + Get_XC_Potential(charge_density.Spin_Summed_Data.vec[i])).ToString());

            // and multiply the density by -e to get the charge density (as these are electrons)
            density = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data density, Band_Data chem_pot)
        {
            // artificially deepen the copies of spin up and spin down
            Band_Data tmp_spinup = new Band_Data(new DoubleMatrix(density.Spin_Up.mat.Rows, density.Spin_Up.mat.Cols));
            Band_Data tmp_spindown = new Band_Data(new DoubleMatrix(density.Spin_Down.mat.Rows, density.Spin_Down.mat.Cols));

            for (int i = 0; i < density.Spin_Up.mat.Rows; i++)
                for (int j = 0; j < density.Spin_Up.mat.Cols; j++)
                {
                    tmp_spinup.mat[i, j] = density.Spin_Up.mat[i, j];
                    tmp_spindown.mat[i, j] = density.Spin_Down.mat[i, j];
                }

            SpinResolved_Data new_density = new SpinResolved_Data(tmp_spinup, tmp_spindown);

            // finally, get the charge density and send it to this new array
            Get_ChargeDensity(layers, ref new_density, chem_pot);

            return new_density;
        }

        public override double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input)
        {
            throw new NotImplementedException();
        }

        public override void Close()
        {
            Console.WriteLine("Closing density solver");
        }

        DoubleHermitianMatrix Create_Hamiltonian(ILayer[] layers, SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(ny * nz);

            // set off diagonal elements 
            for (int i = 0; i < ny - 1; i++)
                for (int j = 0; j < nz - 1; j++)
                {
                    // coupling sites in the growth direction
                    result[i * nz + j + nz, i * nz + j] = ty; result[i * nz + j, i * nz + j + nz] = ty;
                    // coupling sites in the transverse direction
                    result[i * nz + j + 1, i * nz + j] = tz; result[i * nz + j, i * nz + j + 1] = tz;
                }

            double[,] potential = new double[ny, nz];
            // set diagonal elements
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    potential[i, j] = chem_pot.mat[i, j];
                    result[i * nz + j, i * nz + j] = -2.0 * ty + -2.0 * tz + potential[i, j];
                }

            return result;
        }

        double Density_Of_States(double E)
        {
            return Math.Min(Math.Sqrt(Physics_Base.mass / (Physics_Base.hbar * Physics_Base.hbar * E * Math.PI * Math.PI)) * Get_Fermi_Function(E)
                        * DoubleComplex.Norm(tmp_eigvec[tmp_yval * nz + tmp_zval]) * DoubleComplex.Norm(tmp_eigvec[tmp_yval * nz + tmp_zval]), dos_cutoff);
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
    }
}
