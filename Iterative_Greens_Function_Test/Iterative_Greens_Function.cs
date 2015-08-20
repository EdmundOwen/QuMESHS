using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Solver_Bases;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;
using System.IO;

namespace Iterative_Greens_Function_Test
{
    public enum Boundary
    {
        left,
        right
    }

    class Iterative_Greens_Function
    {
        IExperiment exp;

        double dx, dy;
        int nx, ny;

        double alphax, alphay;

        DoubleHermitianMatrix[] onsite_G;
        DoubleHermitianMatrix[] hopping_G;

        /// <summary>
        /// Error for calculating the propagating modes for the boundary conditions.
        /// Their eigenvalue should be unity but this is the window of error
        /// </summary>
        double propagating_mode_error = 1e-3;

        public Iterative_Greens_Function(IExperiment exp)
        {
            this.exp = exp;

          //  this.dx = exp.Dx_Dens; this.dy = exp.Dy_Dens;
          //  this.nx = exp.Nx_Dens; this.ny = exp.Ny_Dens;

            dx = 5.0; dy = 5.0;
            nx = 100; ny = 20;

            alphax = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            alphay = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);







            nx = 100; ny = 1;

            alphax = -1.0;// -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            alphay = 0.0;// -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);
        }

        public void Iterate()
        {
            // cycle through the system at a given energy
            double energy = -0.1;
            Get_Greens_Functions(energy);

        }

        void Get_Greens_Functions(double energy)
        {
            Initiate_Greens_Function(energy);
            for (int i = 1; i < nx - 1; i++)
                Add_Slice(energy, i);

            Output_DoS();
        }

        void Output_DoS()
        {
            StreamWriter sw = new StreamWriter("DoS.dat");
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                    sw.Write((-1.0 * onsite_G[j][i, i].Imag / Math.PI).ToString() + '\t');
                sw.WriteLine();
            }
            sw.Close();
        }

        void Initiate_Greens_Function(double energy)
        {
            onsite_G = new DoubleHermitianMatrix[nx];
            hopping_G = new DoubleHermitianMatrix[nx];

            for (int i = 1; i < nx - 1; i++)
            {
                onsite_G[i] = new DoubleHermitianMatrix(ny);
                hopping_G[i] = new DoubleHermitianMatrix(ny);
            }

            onsite_G[0] = Calculate_Boundary_Conditions(Boundary.left, energy);
            onsite_G[nx - 1] = Calculate_Boundary_Conditions(Boundary.right, energy);
        }

        void Add_Slice(double energy, int slice_no)
        {
            // work out where this slice is
            double x = (slice_no - (nx - 1) / 2) * dx;

            // generate the matrices for this slice
            DoubleHermitianMatrix new_slice = Generate_Slice_Hamiltonian(x);
            DoubleHermitianMatrix new_hopping = Generate_Hopping_Hamiltonian();
            DoubleHermitianMatrix new_hopping_trans = Generate_Hopping_Hamiltonian().Transpose();
            DoubleHermitianMatrix temp_matrix = new DoubleHermitianMatrix(energy * DoubleMatrix.Identity(ny)) - new_slice - Product(Product(new_hopping, onsite_G[slice_no - 1]), new_hopping_trans);

            // factorise the matrix
            DoubleHermitianFact factorisation = new DoubleHermitianFact(temp_matrix);

            // the inverse of this matrix is the onsite Green's function for the new slice
            onsite_G[slice_no] = factorisation.Inverse();

            // the first hop from n to n + 1 uses the onsite Green's function as G_(i,n+1) = G_(i,n) H' G(n+1,n+1) and when i = n, G_(i,n) = G(i,i)
            // must do this here before updating the onsite Green's functions
            hopping_G[slice_no - 1] = onsite_G[slice_no - 1];

            // update the old onsite Green's functions (not including the boundary slice)
            for (int i = 1; i < slice_no; i++)
                onsite_G[i] = onsite_G[i] + Product(Product(Product(Product(hopping_G[i], new_hopping), onsite_G[slice_no]), new_hopping_trans), new DoubleComplexMatrix(MatrixFunctions.ToGeneralMatrix(hopping_G[i])).Transpose());

            // calculate the new hopping Green's functions
            for (int i = 0; i < slice_no; i++)
                hopping_G[i] = Product(Product(hopping_G[i], new_hopping), onsite_G[slice_no]);
        }

        private DoubleHermitianMatrix Calculate_Boundary_Conditions(Boundary boundary, double energy)
        {
            DoubleComplexMatrix bc_matrix = new DoubleComplexMatrix(2 * ny, 2 * ny);

            // generate Hamiltonians for the boundaries
            DoubleHermitianMatrix tmp_slice;
            if (boundary == Boundary.left)
                tmp_slice = Generate_Slice_Hamiltonian(-0.5 * (nx - 1) * dx);
            else if (boundary == Boundary.right)
                tmp_slice = Generate_Slice_Hamiltonian(0.5 * (nx - 1) * dx);
            else throw new NotImplementedException();
            DoubleHermitianMatrix tmp_hopping_trans = Generate_Hopping_Hamiltonian().Transpose();

            // create temporary matrix for top-left of BC matrix
            DoubleHermitianMatrix tmp_matrix = Product(tmp_hopping_trans, (new DoubleHermitianMatrix(energy * DoubleMatrix.Identity(ny)) - tmp_slice));

            // fill with transfer matrix
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < ny; j++)
                {
                    bc_matrix[i, j] = tmp_matrix[i, j];
                    bc_matrix[i + ny, j] = tmp_hopping_trans[i, j];
                    bc_matrix[i, j + ny] = -1.0 * tmp_hopping_trans[i, j];
                }

            // solve eigen-problem for transfer matrix
            DoubleComplexEigDecomp eig_decomp = new DoubleComplexEigDecomp(bc_matrix);
            DoubleComplexEigDecompServer server = new DoubleComplexEigDecompServer();
            eig_decomp = server.Factor(bc_matrix);

            // fill the eigenvalue matrix, excluding unphysical solutions which blow-up in the wire
            DoubleComplexMatrix eig_vals = new DoubleComplexMatrix(ny, ny);
            DoubleComplexMatrix eig_vecs = new DoubleComplexMatrix(ny, ny);
            int count = 0;
            for (int i = 0; i < 2 * ny; i++)
                if (Add_BC_EigenSolution(i, eig_decomp, boundary))
                {
                    // invert eigenvalues if calculating for the right lead
                    if (boundary == Boundary.left)
                        eig_vals[count, count] = eig_decomp.EigenValue(i);
                    else if (boundary == Boundary.right)
                        eig_vals[count, count] = 1.0 / eig_decomp.EigenValue(i);
                    else throw new NotImplementedException();

                    // insert corresponding (normalised) eigenvector into the matrix
                    double norm = 0.0;
                    for (int j = 0; j < ny; j++)
                        norm += DoubleComplex.Norm(eig_decomp.RightEigenVector(i)[j]) * DoubleComplex.Norm(eig_decomp.RightEigenVector(i)[j]);
                    double inv_norm = 1.0 / Math.Sqrt(norm);
                    for (int j = 0; j < ny; j++)
                        eig_vecs[j, count] = inv_norm * eig_decomp.RightEigenVector(i)[j];

                    count++;
                }

            // calculate coefficient matrix A and output Green's function
            DoubleComplexMatrix coefficient_matrix = new DoubleComplexMatrix(ny, ny);
            coefficient_matrix = MatrixFunctions.ToGeneralMatrix(Product(Product(tmp_hopping_trans.Transpose(), eig_vecs), eig_vals));
            DoubleComplexLUFact lu_fact = new DoubleComplexLUFact(coefficient_matrix);
            return new DoubleHermitianMatrix(Product(eig_vecs, lu_fact.Inverse()));

            /*          // allocate for evanescent modes
                      if (boundary == Boundary.left && DoubleComplex.Norm(eig_decomp.EigenValue(i)) < 1.0 - propagating_mode_error)
                          eig_vals[i, i] = eig_decomp.EigenValue(i);
                      else if (boundary == Boundary.right && DoubleComplex.Norm(eig_decomp.EigenValue(i)) > 1.0 + propagating_mode_error)
                          eig_vals[i, i] = 1.0 / eig_decomp.EigenValue(i);
                      // and for propagating modes
                      else if (Math.Abs(DoubleComplex.Norm(eig_decomp.EigenValue(i)) - 1.0) < propagating_mode_error)
                          eig_vals[i, i]
                      else
                          throw new InvalidArgumentException("Error - Cannot have eigenvalue of transfer matrix with value " + eig_decomp.EigenValue(i).ToString() + "!");

                  // fill the eigenvector matrix with only the top half of the eigenvector
                  DoubleComplexMatrix eig_vecs = new DoubleComplexMatrix(ny, ny);
                  DoubleComplexMatrix inv_eigvec = new DoubleComplexMatrix(ny, ny);
                  for (int i = 0; i < ny; i++)
                  {
                      DoubleComplexVector tmpvec_right = eig_decomp.RightEigenVector(i);
                      DoubleComplexVector tmpvec_left = eig_decomp.LeftEigenVector(i);

                      // normalise the top of half of the eigenvector
                      double norm2_right = 0.0;
                      double norm2_left = 0.0;
                      for (int j = 0; j < ny; j++)
                      {
                          norm2_right += NMathFunctions.Abs(tmpvec_right[i]) * NMathFunctions.Abs(tmpvec_right[i]);
                          norm2_left += NMathFunctions.Abs(tmpvec_left[i]) * NMathFunctions.Abs(tmpvec_left[i]);
                      }
                      double norm_left = Math.Sqrt(norm2_left);
                      double norm_right = Math.Sqrt(norm2_right);

                      // and insert it into the matrix
                      for (int j = 0; j < ny; j++)
                      {
                          eig_vecs[j, i] = tmpvec_right[j] / norm_right;
                          inv_eigvec[i, j] = tmpvec_left[j] / norm_left;
                      }
                  }

                  // get the inverse of the eigenvector matrix
       //           DoubleComplexMatrix tmp_eigvec = new DoubleComplexMatrix(eig_vecs);
       //           DoubleComplexMatrix inv_eigvec = NMathFunctions.PseudoInverse(tmp_eigvec);
                  //DoubleComplexLUFact eigvec_fact = new DoubleComplexLUFact(tmp_eigvec);
                  //DoubleComplexMatrix inv_eigvec = eigvec_fact.Inverse();

                  // Calculate the on-site Greens function of the end of the wire and return it
                  if (boundary == Boundary.left)
                      return Product(Product(Product(eig_vecs, eig_vals), inv_eigvec), tmp_hopping_trans);
                  else if (boundary == Boundary.right)
                      return Product(Product(Product(eig_vecs, eig_vals), inv_eigvec), tmp_hopping);
                  else
                      throw new NotImplementedException();
                  */
        }

        /// <summary>
        /// checks to see whether the solution number eigenvalue and eigenvector is unphysical and therefore shouldn't be added
        /// to the boundary condition matrix
        /// </summary>
        bool Add_BC_EigenSolution(int solution_no, DoubleComplexEigDecomp eig_decomp, Boundary boundary)
        {
            double eigenvalue_norm = DoubleComplex.Norm(eig_decomp.EigenValue(solution_no));

            // check for propagating solutions
            if (Math.Abs(eigenvalue_norm - 1.0) < propagating_mode_error)
            {
                throw new NotImplementedException();
            }
            // for evanescent solutions, check for unphysical solutions, ie. |alpha| > 1 for left lead
            else if (eigenvalue_norm > 1.0 && boundary == Boundary.left)
                return true;
            // and |alpha| < 1 for right lead
            else if (eigenvalue_norm < 1.0 && boundary == Boundary.right)
                return true;

            // if these properties are not satisfied, do not add the eigenvalue and eigenvector to the BC (default behaviour)
            return false;
        }

        void Save_BC_Eigvecs(DoubleComplexMatrix eig_mat)
        {
            StreamWriter sw = new StreamWriter("vecs.dat");
            for (int i = 0; i < eig_mat.Rows; i++)
            {
                for (int j = 0; j < eig_mat.Cols; j++)
                    sw.Write(eig_mat[i, j].Real.ToString() + '\t');
                sw.WriteLine();
            }
            sw.Close();
        }

        DoubleHermitianMatrix Generate_Hopping_Hamiltonian()
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(ny);
            for (int i = 0; i < ny; i++)
                result[i, i] = alphax;

            return result;
        }

        DoubleHermitianMatrix Generate_Slice_Hamiltonian(double x)
        {
            DoubleHermitianMatrix result = new DoubleHermitianMatrix(ny);
            for (int i = 0; i < ny - 1; i++)
            {
                double y = (i - (ny - 1) / 2) * dy;

                result[i, i] = -2.0 * alphay - 2.0 * alphax + Get_Potential(x, y);
                result[i, i + 1] = 1.0 * alphay;
                result[i + 1, i] = 1.0 * alphay;
            }
            // the alphax term is necessary here as we are doing this in 2D
            result[ny - 1, ny - 1] = -2.0 * alphax - 2.0 * alphay + Get_Potential(x, (ny - 1) / 2 * dy);

            return result;
        }

        double split_width = 400;
        double split_length = 200;
        double V_max = 0;
        double Get_Potential(double x, double y)
        {
            return V_max * Convert.ToDouble((Math.Abs(x) < 0.5 * split_length) && (Math.Abs(y) > 0.5 * split_width));
        }

        DoubleHermitianMatrix Product(DoubleHermitianMatrix A, DoubleHermitianMatrix B)
        {
            return new DoubleHermitianMatrix(NMathFunctions.Product(MatrixFunctions.ToGeneralMatrix(A), MatrixFunctions.ToGeneralMatrix(B)));
        }

        DoubleHermitianMatrix Product(DoubleComplexMatrix A, DoubleHermitianMatrix B)
        {
            return new DoubleHermitianMatrix(NMathFunctions.Product(A, MatrixFunctions.ToGeneralMatrix(B)));
        }

        DoubleHermitianMatrix Product(DoubleHermitianMatrix A, DoubleComplexMatrix B)
        {
            return new DoubleHermitianMatrix(NMathFunctions.Product(MatrixFunctions.ToGeneralMatrix(A), B));
        }

        DoubleComplexMatrix Product(DoubleComplexMatrix A, DoubleComplexMatrix B)
        {
            return NMathFunctions.Product(A, B);
        }
    }
}
