using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using CenterSpace.NMath.Analysis;

namespace OneD_ThomasFermiPoisson
{
    class OneD_DFTSolver : Density_Base
    {
        double no_kb_T = 20;    // number of kb_T to integrate to
        double t;
        int max_wavefunction = 0;

        int tmp_zval;
        Double tmp_eigval;
        DoubleComplexVector tmp_eigvec;

        public OneD_DFTSolver(double temperature, double dz, int nz, double zmin)
            : base(temperature, 1.0, 1.0, dz, 1, 1, nz, 0.0, 0.0, zmin)
        {
            t = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dz * dz);
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            DoubleHermitianMatrix hamiltonian = Create_Hamiltonian(layers, charge_density, chem_pot);
            DoubleHermitianEigDecomp eig_decomp = new DoubleHermitianEigDecomp(hamiltonian);

            double min_eigval = eig_decomp.EigenValues.Min();
            max_wavefunction = (from val in eig_decomp.EigenValues
                                    where val < no_kb_T * Physics_Base.kB * temperature
                                    select val).ToArray().Length;

            DoubleVector dens_up = new DoubleVector(nz, 0.0);
            DoubleVector dens_down = new DoubleVector(nz, 0.0);
            OneVariableFunction dens_of_states = new OneVariableFunction(new Func<double, double>(Density_Of_States));

            for (int j = 0; j < nz; j++)
            {
                double dens_val = 0.0;
                for (int i = 0; i < max_wavefunction; i++)
                {
                    // set temporary eigenvalue and eigenvector
                    tmp_eigval = eig_decomp.EigenValue(i); tmp_eigvec = eig_decomp.EigenVector(i);
                    // and position
                    tmp_zval = j;

                    // and integrate the density of states at this position for this eigenvector from the minimum energy to
                    // (by default) 20 * k_b * T above mu = 0
                    dens_val += dens_of_states.Integrate(min_eigval, no_kb_T * Physics_Base.kB * temperature);
                }

                // just share the densities (there is no spin polarisation)
                dens_up[j] = 0.5 * dens_val;
                dens_down[j] = 0.5 * dens_val;
            }

            //
            //Console.Clear();
            //for (int i = 0; i < nz; i++)
            //    Console.WriteLine(charge_density.Spin_Summed_Data[i].ToString() + "\t" + (-1.0 * Physics_Base.q_e * (dens_up[i] + dens_down[i])).ToString() + "\t" + (chem_pot.vec[i] + Get_XC_Potential(charge_density.Spin_Summed_Data.vec[i])).ToString());

            // and multiply the density by -e to get the charge density (as these are electrons)
            charge_density = -1.0 * Physics_Base.q_e * new SpinResolved_Data(new Band_Data(dens_up), new Band_Data(dens_down));
        }

        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data density, Band_Data chem_pot)
        {
            SpinResolved_Data new_density = new SpinResolved_Data(new Band_Data(density.Spin_Up.vec.DeepenThisCopy()), new Band_Data(density.Spin_Down.vec.DeepenThisCopy()));

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
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(nz);

            // set off diagonal elements
            for (int i = 0; i < nz - 1; i++)
            {
                result[i + 1, i] = t; result[i, i + 1] = t;
            }

            double[] potential = new double[nz];
            // set diagonal elements
            for (int i = 0; i < nz; i++)
            {
                potential[i] = chem_pot.vec[i];
                result[i, i] = -2.0 * t + potential[i];
            }

            return result;
        }

        double Density_Of_States(double E)
        {
            return Physics_Base.mass / (Physics_Base.hbar * Physics_Base.hbar * 2.0 * Math.PI) * Get_Fermi_Function(E)
                        * DoubleComplex.Norm(tmp_eigvec[tmp_zval]) * DoubleComplex.Norm(tmp_eigvec[tmp_zval]);
        }

        public int No_Wavefunctions
        {
            get { return max_wavefunction; }
        }
    }
}
