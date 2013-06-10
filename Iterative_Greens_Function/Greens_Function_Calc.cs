using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using OpenCLNet;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;

namespace Iterative_Greens_Function
{
    class Greens_Function_Calc
    {
        int no_Slices, slice_width;
        double dx;
        DoubleComplexMatrix hop_mat, hop_mat_conj, H0_slice;

        // physical constants
        public const double hbar = 0.658211814;             // (meV) (ps)
        public const double q_e = -0.160217646;             // (aC)
        public const double mass = 0.067 * 5.68562958e-3;    // (meV) (ps)^2 (nm)^-2 with GaAs effective mass
        public const double epsilon = 13 * 1.41859713e-6;   // (aC)^2 (nm)^-1 (meV)^-1 for GaAs

        // dimensionality parameter
        double alpha;

        public Greens_Function_Calc(double dx, int no_Slices, int slice_width)
        {
            this.dx = dx;
            this.no_Slices = no_Slices;
            this.slice_width = slice_width;

            // Hamiltonian hopping prefactor  
            // (all energies are divided by this factor)
            alpha = hbar * hbar / (2.0 * mass * dx * dx);
        }

        /// <summary>
        /// Initialise Green's function iterator
        /// </summary>
        public void Initialise()
        {
            Console.WriteLine("Initialising Green's function iterator");

            // create hopping matrices (dimensionless)
            hop_mat = new DoubleComplexMatrix(slice_width, slice_width);
            for (int i = 0; i < slice_width; i++)
                    hop_mat[i, i] = -1.0;
            hop_mat_conj = NMathFunctions.ConjTranspose(hop_mat);

            // create on-slice Hamiltonian without potential (dimensionless)
            H0_slice = new DoubleComplexMatrix(slice_width, slice_width);
            for (int i = 0; i < slice_width; i++)
            {
                H0_slice[i, i] = 2.0;

                if (i != 0)
                    H0_slice[i, i - 1] = -1.0;
                if (i != slice_width - 1)
                    H0_slice[i, i + 1] = -1.0;
            }
        }

        public DoubleComplexMatrix[] Iterate(out DoubleComplexMatrix G_0n, out DoubleComplexMatrix G_n0, DoubleMatrix Potential, double Energy)
        {
            DoubleComplexMatrix[] G_ii = new DoubleComplexMatrix[no_Slices];
            DoubleComplexMatrix[] G_in = new DoubleComplexMatrix[no_Slices];
            DoubleComplexMatrix[] G_ni = new DoubleComplexMatrix[no_Slices];

            // normalise energy and potential
            Energy = Energy / alpha;
            Potential = Potential / alpha;

            Console.WriteLine("Calculating initial boundary conditions");
            // Calculate Green's function boundary conditions
            Calculate_BCs(ref G_ii, Potential,  Energy);
            // and connecting Green's function initial conditions are trivial as i = 0, n = 0
            G_in[0] = (DoubleComplexMatrix)G_ii[0].Clone();
            G_ni[0] = (DoubleComplexMatrix)G_ii[0].Clone();

            Console.WriteLine("Iterating slices");
            // Iterate Green's function at given energy
            for (int j = 1; j < no_Slices - 1; j++)
                Iterate_Slice(ref G_ii, ref G_in, ref G_ni, Potential, j, Energy);

            // Add output lead (the Green's function for this slice was calculated
            // with the boundary conditions)
            Calculate_Connections(G_ii[no_Slices - 1], ref G_ii, ref G_in, ref G_ni, Potential, no_Slices - 1, Energy);

            // copy G_0n and G_n0 to the signiture inputs for current calculations
            G_0n = G_in[0];
            G_n0 = G_ni[0];

            return G_ii;
        }

        void Calculate_BCs(ref DoubleComplexMatrix[] G_ii, DoubleMatrix Potential, double Energy)
        {
            // Generate transfer matrix
            DoubleComplexMatrix transfer_mat = Generate_Transfer_Matrix(Potential, Energy);

            Console.WriteLine("Calculating eigen-decomposition of BC transfer matrix");
            // Calculate eigenvalues
            DoubleComplexEigDecomp eig_decomp = new DoubleComplexEigDecomp(transfer_mat);
            if (!eig_decomp.IsGood)
                throw new Exception("Error - Eigenvalue decomposition of transfer matrix failed!");

            Console.WriteLine("Sorting eigenvalues");
            // Sort through eigenvalues
            int[] backward_modes, forward_modes;
            int no_forwards, no_backwards, no_propagating_modes;
            Sort_Eigenvalues(eig_decomp, out backward_modes, out forward_modes, out no_forwards, out no_backwards, out no_propagating_modes);

            Console.WriteLine("Compiling eigenvectors");
            // Compile eigenvectors
            DoubleComplexMatrix forward_eigenvecs, backward_eigenvecs;
            DoubleComplexMatrix forward_eigenvals, backward_eigenvals;
            Compile_Eigenvectors(out forward_eigenvecs, out backward_eigenvecs, out forward_eigenvals, out backward_eigenvals, eig_decomp, backward_modes, forward_modes, no_forwards, no_backwards);
            Normalise_Eigenvectors(ref forward_eigenvecs, ref backward_eigenvecs);

            Console.WriteLine("Creating boundary conditions");
            // Input boundary conditions into G_ii
            G_ii = new DoubleComplexMatrix[no_Slices];
            Create_Boundary_Conditions(ref G_ii, forward_eigenvecs, backward_eigenvecs, forward_eigenvals, backward_eigenvals);
        }

        /// <summary>
        /// Generates the transfer matrix needed to solve the boundary conditions from the leads
        /// NOTE!!!!!! This is calculated for the first slice of the potential which assumes that
        /// the system is symmetric!!
        /// </summary>
        DoubleComplexMatrix Generate_Transfer_Matrix(DoubleMatrix Potential, double Energy)
        {
            DoubleComplexMatrix trans = new DoubleComplexMatrix(2 * slice_width, 2 * slice_width, 0);

            ////////////// NOTE: Hamiltonian is calculated for the potential on the first slice!!!!!! ////////////
            DoubleComplexMatrix top_left = NMathFunctions.Product(hop_mat_conj, Calculate_Hamiltonian(Potential, 0, Energy));
            for (int i = 0; i < slice_width; i++)
                for (int j = 0; j < slice_width; j++)
                {
                    // fill top left
                    trans[i, j] = top_left[i, j];

                    // fill off diagonal blocks
                    trans[i + slice_width, j] = hop_mat_conj[i, j];
                    trans[i, j + slice_width] = -1 * hop_mat_conj[i, j];
                }

            return trans;
        }

        /// <summary>
        /// sorts the eigenvalues to work out which have |alpha| \geq 1 and which have
        /// |alpha| \leq 1 to discount unphysical solutions in the leads
        /// </summary>
        void Sort_Eigenvalues(DoubleComplexEigDecomp eig_decomp, out int[] backward_modes, out int[] forward_modes, out int no_forward, out int no_backward, out int no_propagating)
        {
            // a tolerance for the propagating eigenvalues around |alpha| = 1
            double tol = 0.00000001;

            no_forward = 0;
            no_backward = 0;
            no_propagating = 0;

            DoubleVector abs_vals = NMathFunctions.Abs(eig_decomp.EigenValues);
            forward_modes = new int[abs_vals.Length];
            backward_modes = new int[abs_vals.Length];

            // and add the assign
            for (int i = 0; i < abs_vals.Length; i++)
            {
                // first check for propagating modes
                if (abs_vals[i] >= 1.0 - tol && abs_vals[i] <= 1.0 + tol)
                {
                    double current = Check_Direction_Of_Propagation(eig_decomp, i);
                    if (current > 0)
                    {
                        forward_modes[no_forward] = i;
                        no_forward++;
                        no_propagating++;
                    }
                    else if (current < 0)
                    {
                        backward_modes[no_backward] = i;
                        no_backward++;
                        no_propagating++;
                    }
                    else
                        throw new Exception("Error - Propagating mode is not propagating!");
                }

                // then add the evanescent modes
                if (abs_vals[i] < 1.0 - tol)
                {
                    forward_modes[no_forward] = i;
                    no_forward++;
                }
                else if (abs_vals[i] > 1.0 + tol)
                {
                    backward_modes[no_backward] = i;
                    no_backward++;
                }
            }
        }

        /// <summary>
        /// Calculates the direction of propagation of the i-th mode from the matrix decomposition
        /// Taken from Crispin's code.... ie. Check how this works!
        /// </summary>
        double Check_Direction_Of_Propagation(DoubleComplexEigDecomp eig_decomp, int i)
        {
            double current = 0.0;

            DoubleComplex eigenvalue = eig_decomp.EigenValue(i);
            DoubleComplexVector eigenvector = eig_decomp.RightEigenVector(i);
            for (int j = 0; j < slice_width; j++)
            {
                double tmp = NMathFunctions.Abs(eigenvector[j]);
                current += 2.0 * (eigenvalue * hop_mat[j, j]).Imag * tmp * tmp;
            }

            return current;
        }

        /// <summary>
        /// Compiles the eigenvectors in either direction from the eigenvalue decomposition
        /// using the results from "Sort_Eigenvalues"
        /// </summary>
        void Compile_Eigenvectors(out DoubleComplexMatrix forward_eigenvecs, out DoubleComplexMatrix backward_eigenvecs, out DoubleComplexMatrix forward_eigenvals, out DoubleComplexMatrix backward_eigenvals, 
            DoubleComplexEigDecomp eig_decomp, int[] backward_modes, int[] forward_modes, int no_forwards, int no_backwards)
        {
            forward_eigenvecs = new DoubleComplexMatrix(slice_width, slice_width, 0.0);
            forward_eigenvals = new DoubleComplexMatrix(slice_width, slice_width, 0.0);
            for (int i = 0; i < no_forwards; i++)
            {
                forward_eigenvals[i, i] = eig_decomp.EigenValue(i);
                for (int j = 0; j < slice_width; j++)
                    forward_eigenvecs[j, i] = eig_decomp.RightEigenVector(forward_modes[i])[j];
            }

            backward_eigenvecs = new DoubleComplexMatrix(slice_width, slice_width, 0.0);
            backward_eigenvals = new DoubleComplexMatrix(slice_width, slice_width, 0.0);
            for (int i = 0; i < no_backwards; i++)
            {
                backward_eigenvals[i, i] = eig_decomp.EigenValue(i);
                for (int j = 0; j < slice_width; j++)
                    backward_eigenvecs[j, i] = eig_decomp.RightEigenVector(backward_modes[i])[j];
            }
        }

        /// <summary>
        /// Normalises the eigenvector matrices
        /// </summary>
        void Normalise_Eigenvectors(ref DoubleComplexMatrix forward_eigenvecs, ref DoubleComplexMatrix backward_eigenvecs)
        {
            // normalise each column of the forward eigenvector matrix
            for (int i = 0; i < forward_eigenvecs.Cols; i++)
            {
                double norm2 = forward_eigenvecs.Slice(0, i, forward_eigenvecs.Rows, 1, 0).TwoNorm();
                double inv_norm = 1.0 / norm2;

                for (int j = 0; j < forward_eigenvecs.Rows; j++)
                    forward_eigenvecs[j, i] = forward_eigenvecs[j, i] * inv_norm;
            }

            // and similarly for backwards eigenvectors
            for (int i = 0; i < backward_eigenvecs.Cols; i++)
            {
                double norm2 = backward_eigenvecs.Slice(0, i, backward_eigenvecs.Rows, 1, 0).TwoNorm();
                double inv_norm = 1.0 / norm2;

                for (int j = 0; j < backward_eigenvecs.Rows; j++)
                    backward_eigenvecs[j, i] = backward_eigenvecs[j, i] * inv_norm;
            }
        }

        /// <summary>
        /// Calculates and inputs the boundary conditions into the first and last elements of G_ii
        /// </summary>
        void Create_Boundary_Conditions(ref DoubleComplexMatrix[] G_ii, DoubleComplexMatrix forward_eigenvecs, DoubleComplexMatrix backward_eigenvecs, DoubleComplexMatrix forward_eigenvals, DoubleComplexMatrix backward_eigenvals)
        {
            G_ii[0] = NMathFunctions.Product(NMathFunctions.Product(NMathFunctions.Product(forward_eigenvecs, forward_eigenvals), NMathFunctions.Inverse(forward_eigenvecs)), hop_mat_conj);

            Print_Matrix("forward_eigenvecs", forward_eigenvecs);
            Print_Matrix("forward_eigenvals", forward_eigenvals);
            Print_Matrix("G_ii", G_ii[0]);
            Print_Matrix("backward_eigenvals", backward_eigenvals);
            Print_Matrix("backward_eigenvecs", backward_eigenvecs);

            G_ii[no_Slices - 1] = NMathFunctions.Product(NMathFunctions.Product(NMathFunctions.Product(backward_eigenvecs, NMathFunctions.Inverse(backward_eigenvals)), NMathFunctions.Inverse(backward_eigenvecs)), hop_mat);
        }
        
        void Iterate_Slice(ref DoubleComplexMatrix[] G_ii, ref DoubleComplexMatrix[] G_in, ref DoubleComplexMatrix[] G_ni, DoubleMatrix Potential, int slice_no, double Energy)
        {
            // Calculate G_n+1,n+1 (MacKinnon 5a)
            DoubleComplexMatrix E_minus_H = Calculate_Hamiltonian(Potential, slice_no, Energy);
            G_ii[slice_no] = Calculate_New_Greens_Function(E_minus_H, G_ii[slice_no - 1]);
            
            // Calculate influence of the new slice on the previously calculated results
            Calculate_Connections(G_ii[slice_no], ref G_ii, ref G_in, ref G_ni, Potential, slice_no, Energy);
        }
        
        /// <summary>
        /// Calculates the Green's functions from each slice to the end and 
        /// how the new diagonal Green's functions are influenced by the new
        /// slice
        /// </summary>
        void Calculate_Connections(DoubleComplexMatrix G_nplusnplus, ref DoubleComplexMatrix[] G_ii, ref DoubleComplexMatrix[] G_in, ref DoubleComplexMatrix[] G_ni, DoubleMatrix Potential, int slice_no, double Energy)
        {
            if (slice_no % 10 == 0)
                Console.WriteLine("Calculating slice no:\t" + slice_no.ToString());

            // Calculate G_n+1,i 
            G_ni[slice_no] = (DoubleComplexMatrix)G_nplusnplus.Clone();
            for (int i = 0; i < slice_no; i++)
                G_ni[i] = NMathFunctions.Product(NMathFunctions.Product(G_ii[slice_no], hop_mat_conj), G_ni[i]);// ViennaCL.TripMatProd(G_ii[slice_no], hop_mat_conj, G_ni[i]);   // (MacKinnon 5c)

            // Calculate G_i,i ( note that G_ni has been updated to include the next slice so
            // it is technically G_{n+1, i} )
            for (int i = 1; i < slice_no; i++)
                G_ii[i] = G_ii[i] + NMathFunctions.Product(NMathFunctions.Product(G_in[i], hop_mat), G_ni[i]); // ViennalCL.Add(G_ii[i], ViennalCL.TripProd(G_in[i], hop_mat, G_ni[i]));      // (MacKinnon 5d with i = j)

            // Calculate G_i,n+1
            G_in[slice_no] = (DoubleComplexMatrix)G_nplusnplus.Clone();
            for (int i = 1; i < slice_no + 1; i++)
                G_in[i] = NMathFunctions.Product(NMathFunctions.Product(G_in[i], hop_mat), G_ii[slice_no]); // ViennaCL.TripProd(G_in[i], hop_mat, G_ii[slice_no]);  // (MacKinnon 5d)
        }

        DoubleComplexMatrix Calculate_Hamiltonian(DoubleMatrix Potential, int slice_no, double Energy)
        {
            DoubleComplexMatrix result = new DoubleComplexMatrix(slice_width, slice_width);
            result.Diagonal().Set(Range.All, Energy);

            // Get potential from memory
            DoubleVector pot_data = Potential.Slice(0, slice_no, slice_width, 1, 0);

            DoubleMatrix pot_slice = new DoubleMatrix(slice_width, slice_width);// ViennaCL.GetSlice(Potential, slice_start, slice_start + slice_width, 1);
            for (int i = 0; i < slice_width; i++)
                pot_slice[i, i] = Potential[i, slice_no];

            Print_Matrix("pot_slice", pot_slice);
            
            // Construct (E - V) - H_0
            return result - (H0_slice + pot_slice);
        }

        DoubleComplexMatrix Calculate_New_Greens_Function(DoubleComplexMatrix E_minus_H, DoubleComplexMatrix G_nn)
        {
            // Calculate V G_nn V\dag
            DoubleComplexMatrix result = NMathFunctions.Product(NMathFunctions.Product(hop_mat, G_nn), hop_mat_conj);
            result = E_minus_H + result;

            // Invert matrix
            return MatrixFunctions.Inverse(result);
        }

        public double Calculate_Currents()
        {
            throw new NotImplementedException();
        }

        public void Print_Matrix(string filename, DoubleMatrix mat)
        {
            Print_Matrix(filename, (DoubleComplexMatrix)mat);
        }

        void Print_Matrix(string filename, DoubleComplexMatrix mat)
        {
            StreamWriter sw_re = new StreamWriter(filename + "_re");
            StreamWriter sw_im = new StreamWriter(filename + "_im");

            for (int i = 0; i < mat.Rows; i++)
                for (int j = 0; j < mat.Cols; j++)
                {
                    sw_re.Write(mat[i, j].Real + "\t");
                    sw_im.Write(mat[i, j].Imag + "\t");
                    if (j == mat.Cols - 1)
                    {
                        sw_re.WriteLine();
                        sw_im.WriteLine();
                    }
                }

            sw_re.Close();
            sw_im.Close();
        }
    }
}
