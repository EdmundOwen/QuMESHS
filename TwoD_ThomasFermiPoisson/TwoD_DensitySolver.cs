using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    class TwoD_DensitySolver : Density_Solver
    {
        DoubleComplexMatrix H;

        public TwoD_DensitySolver(double dx, double fermi_Energy, int nx) 
            : base(fermi_Energy, 0.0, dx, 1.0, 1.0, nx, 1, 1)
        {
            this.dx = dx; this.nx = nx;
        }

        /// <summary>
        /// Initialises the Hamiltonian with the given potential from the well
        /// </summary>
        public void Initialise_Hamiltonian(DoubleVector well_potential)
        {
            // matrix is 2 * nx due to spin
            H = new DoubleComplexMatrix(2 * nx, 2 * nx, 0.0);

            // this is the hopping element
            double alpha = hbar * hbar / (2.0 * mass * dx * dx);

            // initialise matrix elements
            for (int i = 0; i < 2 * nx - 1; i++)
            {
                // off-diagonal
                H[i + 1, i] = - 1.0 * alpha;
                H[i, i + 1] = - 1.0 * alpha;

                // diagonal modulo number of lattice sites for spin
                H[i, i] = 2.0 * alpha + well_potential[i % nx];
            }
            // and fill the final element
            H[2 * nx - 1, 2 * nx - 1] = 2.0 * alpha + well_potential[nx - 1];

            // and delete the middle elements to make this block diagonal
            // this decouples the nx-1 spin-up component and the 0 spin-down component
            H[nx - 1, nx] = 0.0;
            H[nx, nx - 1] = 0.0;
        }

        public DoubleVector Solve_Density_Using_Momentum(double fermi_Energy, double no_kB_T, double dk)
        {
            // calculate the maximum energy to calculate to above the fermi surface
            double max_Energy = fermi_Energy + (no_kB_T * kB * temperature);

            // density vector is spin resolved
            DoubleVector density = new DoubleVector(2 * nx);
            
            // number of subbands with energy still below the energy we're calculating up to
            int no_subbands = 0;
            // number of k points calculated for
            int no_ky = 0;

            // calculate states and energies at k = 0
            DoubleComplexMatrix k_mat = new DoubleComplexMatrix(2 * nx, 2 * nx);
            k_mat.Diagonal().Set(Range.All, 0.0);
            DoubleComplexEigDecomp eigs = new DoubleComplexEigDecomp();
            eigs.FactorNoPreconditioning(H + k_mat);

            DoubleComplexVector energies_before = eigs.EigenValues;
            DoubleComplexMatrix states_before = eigs.RightEigenVectors;
            Sort_Eigenvectors(ref energies_before, ref states_before);

            // calculate number of energy levels below max energy
            for (int i = 0; i < 2 * nx; i++)
                if (energies_before[i].Real < max_Energy)
                    no_subbands++;

            double cutoff = 0.000001;
            Console.WriteLine("DEBUG: CUTOFF FOR DENSITY OF STATES IMPLEMENTED = " + cutoff.ToString());

            // cycle through momentum until fermi level is reached for all subbands
            while (no_subbands != 0)
            {
                double ky = no_ky * dk;

                // diagonalise H + k_y^2
                k_mat.Diagonal().Set(Range.All, ky * ky);
                eigs.FactorNoPreconditioning(H + k_mat);

                // create new energies and eigenvalues for this value of ky
                DoubleComplexVector new_eig_vals = eigs.EigenValues;
                DoubleComplexMatrix new_eig_states = eigs.RightEigenVectors;
                // and sort them
                Sort_Eigenvectors(ref new_eig_vals, ref new_eig_states);

                //DEBUGGING LINE
                DoubleComplexVector test = new_eig_vals - energies_before;

                // calculate dk/dE for remaining subbands
                DoubleVector[] dn_dk = new DoubleVector[no_subbands];
                double[] dk_dE = new double[no_subbands];
                double[] fermi_function = new double[no_subbands];

                int tmp = no_subbands;
                for (int i = 0; i < tmp; i++)
                    // check that the indexed band is still below the maximum energy
                    if (new_eig_vals[i].Real > max_Energy)
                        no_subbands--;
                    // and if it is, calculate dk/dE and dn(x)/dk and integrate
                    else
                    {
                        // calculate dk/dE (applying a cut off when the bands are flat)
                        double dE = (new_eig_vals[i] - energies_before[i]).Real;
                        if (dE < cutoff)
                            dk_dE[i] = 1.0 / cutoff;
                        else
                            dk_dE[i] = dk / (new_eig_vals[i] - energies_before[i]).Real;

                        // and dn(x) / dk (first term is from the longitudinal direction and second term is the transverse density change, calculated numerically)
                        dn_dk[i] = (2 / Math.PI) * Calculate_Density_Change(states_before, new_eig_states, i);

                        // and get the fermi function for this energy
                        fermi_function[i] = Get_Fermi_Function(new_eig_vals[i].Real);
                    }

                Console.WriteLine((1.0/dk_dE[0]).ToString() + '\t' + (1.0/dk_dE[1]).ToString());

                // Finally, integrate the sub-band densities and add
                for (int i = 0; i < no_subbands; i++)
                    for (int j = 0; j < 2 * nx; j++)
                        density[j] += dn_dk[i][j] * dk_dE[i] * fermi_function[i];

                // and replace the old eigenvectors and eigenvalues
                energies_before = new_eig_vals;
                states_before = new_eig_states;

                // and increment the ky value
                no_ky++;
            }

            return density;
        }

        /// <summary>
        /// Sorts the eigenvectors in ascending energy. Both inputs are changed
        /// </summary>
        void Sort_Eigenvectors(ref DoubleComplexVector energies, ref DoubleComplexMatrix states)
        {
            // First, energies should be real so we can extract this data to an array of doubles
            DoubleVector data = new DoubleVector(energies.ToRealDataTable());
            double[] array_data = data.ToArray();

            int no_elements = array_data.Length;

            // we fill a permutation array with a set of indices
            int[] permutation = new int[no_elements];
            for (int i = 0; i < no_elements; i++)
                permutation[i] = i;

            // and sort them so that the energies are increasing 
            Array.Sort(array_data, permutation);
            // and permute for the spin sectors (see method)
            Permute_Spin_Sectors(ref array_data, ref permutation);

            // we initialise a new vector and matrix for the sorted values
            DoubleComplexVector tmp_energies = new DoubleComplexVector(no_elements);
            DoubleComplexMatrix tmp_states = new DoubleComplexMatrix(no_elements, no_elements);
            // and permute accordingly
            for (int i = 0; i < no_elements; i++)
            {
                // change the energies
                tmp_energies[i] = energies[permutation[i]];

                // and the states
                for (int j = 0; j < no_elements; j++)
                    tmp_states[j, i] = states[j, permutation[i]];
            }

            // and this is the new "energies" vector and "states" matrix
            states = tmp_states;
            energies = tmp_energies;
        }

        /// <summary>
        /// This slightly alters the permutation matrix in the case of spin degeneracy
        /// As the spins are degenerate, we must have a consistent sorting so if the energies are
        /// the same to a given tolerance, the spin-up contribution, $i \leq n_x$ is placed first
        /// </summary>
        void Permute_Spin_Sectors(ref double[] energies, ref int[] permutation)
        {
            double tol = 0.00000001;

            for (int i = 0; i < 2 * nx - 1; i++)
                if (Math.Abs(energies[i] - energies[i + 1]) < tol && permutation[i] > permutation[i + 1])
                {
                    int tmp1 = permutation[i];
                    int tmp2 = permutation[i + 1];
                    // swap these two values
                    permutation[i + 1] = tmp1;
                    permutation[i] = tmp2;

                    // and although the energies should be the same, swap them for consistency
                    double tmp3 = energies[i]; double tmp4 = energies[i + 1];
                    energies[i + 1] = tmp3;
                    energies[i] = tmp4;
                }
        }

        /// <summary>
        /// Calculates the change in density between the states given by the columns of the input matrices with index "sub_band"
        /// </summary>
        private DoubleVector Calculate_Density_Change(DoubleComplexMatrix states_before, DoubleComplexMatrix states_after, int sub_band)
        {
            // Initialise the density vector
            DoubleVector dens_before = new DoubleVector(2 * nx);
            DoubleVector dens_after = new DoubleVector(2 * nx);

            // Calculate the density for the given sub-band (these should be the vertical columns)
            for (int i = 0; i < 2 * nx; i++)
            {
                dens_before[i] = (states_before[i, sub_band] * NMathFunctions.Conj(states_before[i, sub_band])).Real;
                dens_after[i] = (states_after[i, sub_band] * NMathFunctions.Conj(states_after[i, sub_band])).Real;
            }

            DoubleVector result = dens_after - dens_before;

            return result;
        }

        public DoubleVector Solve_Density_Using_Energy(double dE, double init_energy, int Energy_steps)
        {
            // density vector is spin resolved
            DoubleVector density = new DoubleVector(2 * nx);

            // cycle through energy steps up fermi level
            for (int i = 0; i < Energy_steps; i++)
            {
                double Energy = init_energy + i * dE;

                // diagonalise H - E
                DoubleComplexMatrix E_mat = new DoubleComplexMatrix(2 * nx, 2 * nx);
                E_mat.Diagonal().Set(Range.All, Energy);
                DoubleComplexEigDecomp eigs = new DoubleComplexEigDecomp(H - E_mat);

                // calculate density for given energy slice
            }

            throw new NotImplementedException();
            return density;
        }

        public DoubleVector Solve_Density_Using_GreensFunctions(double dE, double init_energy, int Energy_steps)
        {
            // density vector is spin resolved
            DoubleVector density = new DoubleVector(2 * nx);

            // do something

            throw new NotImplementedException();
            return density;
        }
    }
}
