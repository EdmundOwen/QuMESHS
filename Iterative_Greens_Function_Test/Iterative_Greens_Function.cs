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

    public class Iterative_Greens_Function : Density_Base
    {
        IExperiment exp;

        double dx, dy;
        int nx, ny;

        double alphax, alphay, alphax_prime, alphay_prime;

        DoubleComplexMatrix[] onsite_G;
        DoubleComplexMatrix[] hopping_G;

        double[,] potential_xy;

        double temperature;

        DoubleMatrix DensityoS;
        DoubleComplexMatrix LHSboundary;
        DoubleComplexMatrix RHSboundary;


        /// <summary>
        /// Error for calculating the propagating modes for the boundary conditions.
        /// Their eigenvalue should be unity but this is the window of error
        /// </summary>
        double propagating_mode_error = 1e-3;

        public Iterative_Greens_Function(IExperiment exp, double[,] potential)
            : base(exp.Temperature)
        {
            this.exp = exp;

          //  this.dx = exp.Dx_Dens; this.dy = exp.Dy_Dens;
          //  this.nx = exp.Nx_Dens; this.ny = exp.Ny_Dens;

            dx = 100; dy = 100;
            nx = 50; ny = 20;

          alphax_prime = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dx * dx);
          alphay_prime = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dy * dy);

            potential_xy = potential;
            DensityoS = new DoubleMatrix(nx, ny);






           // nx = 100; ny = 1;

            // alphax = -1.0;// -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            // alphay = 0.0;// -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);

            alphax = alphax_prime / alphax_prime;
            alphay = alphay_prime / alphax_prime;
        }

        public Iterative_Greens_Function(IExperiment exp)
            :base (exp.Temperature)
        {
            this.exp = exp;

            //  this.dx = exp.Dx_Dens; this.dy = exp.Dy_Dens;
            //  this.nx = exp.Nx_Dens; this.ny = exp.Ny_Dens;

            dx = 100; dy = 100;
            nx = 50; ny = 20;

            alphax_prime = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dx * dx);
            alphay_prime = -0.5 * Physics_Base.hbar * Physics_Base.hbar / (mass * dy * dy);


            temperature = 4;

            DensityoS = new DoubleMatrix(nx, ny);



            // nx = 100; ny = 1;

            // alphax = -1.0;// -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dx * dx);
            // alphay = 0.0;// -0.5 * Physics_Base.hbar * Physics_Base.hbar / (Physics_Base.mass * dy * dy);

            alphax = alphax_prime / alphax_prime;
            alphay = alphay_prime / alphax_prime;
        }

        public void Iterate()
        {
            // cycle through the system at a given energy
            
            int n_steps = 600;
            double deltaE = 0.001;
            //double emin = -1.0 * deltaE * n_steps / 2;
            double emin = 1.2;
            DoubleVector DoSVector;
            DoubleVector EnergyVector;
            //DoubleVector MaxVector;
            DoSVector = new DoubleVector(n_steps);
            EnergyVector = new DoubleVector(n_steps);
           // MaxVector = new DoubleVector(n_steps);
            
            for (int i = 0; i < n_steps; i++)
            {
                double energy = emin + i * deltaE;
                double epsilon = energy / alphax_prime;
                Get_Greens_Functions(epsilon);
                DoSVector[i] = GetTotalDoS();
                EnergyVector[i] = energy;
               // MaxVector[i] = GetMaxDoS();
            }

            Output_DoSVector(DoSVector, EnergyVector);
            
            
            /*
            double energy = 0.431;
            double epsilon = energy / alphax_prime;
            Get_Greens_Functions(epsilon);

            //bool LHS = (LHSboundary == onsite_G[0]);
            //bool RHS = (RHSboundary == onsite_G[nx - 1]);

            Output_DoS();
            */
        }

        void Get_Greens_Functions(double energy)
        {
            Initiate_Greens_Function(energy);
            for (int i = 1; i < nx - 1; i++)
                Add_Slice(energy, i);
            // factorise and invert the right hand lead Green's function
            DoubleComplexLUFact factorisation = new DoubleComplexLUFact(onsite_G[nx - 1]);
            DoubleComplexMatrix rhs_inverted_G = factorisation.Inverse();
            Bind_Slice(nx - 1, rhs_inverted_G);

            //Output_DoS();
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                    DensityoS[j, i] = (-1.0 * onsite_G[j][i, i].Imag / (Math.PI * alphax_prime));
            }
        }

        double GetTotalDoS()
        {
            double TotalDoS = DensityoS.Average() * nx * ny * dx * dy;

            return TotalDoS;
        }

        double GetMaxDoS()
        {
            double MinDoS = DensityoS.Min();
            double MaxDoS = DensityoS.Max();

            if (MaxDoS > -1.0 * MinDoS)
                return MaxDoS;
            else
                return MinDoS;

        }

        public double[,] GetDoS(double energy)
        {
            double epsilon = energy / alphax_prime;

            Get_Greens_Functions(epsilon);

            double[,] DoS = new double[nx, ny];
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                    DoS[j,i]=(-1.0 * onsite_G[j][i, i].Imag / (Math.PI * alphax_prime));
            }
            return DoS;
        }

        void Output_DoSVector(DoubleVector DoSVector, DoubleVector MaxVector, DoubleVector EnergyVector)
        {
            StreamWriter sw = new StreamWriter("DoSVector.dat");
            for (int i = 0; i < DoSVector.Length; i++)
            {
                sw.Write(EnergyVector[i].ToString() + '\t' + DoSVector[i].ToString() + '\t' + (DoSVector[i]*Get_Fermi_Function(EnergyVector[i], 0.0, temperature)).ToString() + '\t' + MaxVector[i].ToString());
                sw.WriteLine();
            }
            sw.Close();

        }

        void Output_DoSVector(DoubleVector DoSVector, DoubleVector EnergyVector)
        {
            StreamWriter sw = new StreamWriter("DoSVector.dat");
            for (int i = 0; i < DoSVector.Length; i++)
            {
                sw.Write(EnergyVector[i].ToString() + '\t' + DoSVector[i].ToString());
                sw.WriteLine();
            }
            sw.Close();

        }
        void Output_DoS()
        {
            StreamWriter sw = new StreamWriter("DoS.dat");
            for (int i = 0; i < ny; i++)
            {
                for (int j = 0; j < nx; j++)
                    sw.Write((-1.0 * onsite_G[j][i, i].Imag / (Math.PI*alphax_prime)).ToString() + '\t');
                sw.WriteLine();
            }
            sw.Close();
        }

        void Initiate_Greens_Function(double energy)
        {
            onsite_G = new DoubleComplexMatrix[nx];
            hopping_G = new DoubleComplexMatrix[nx];

            for (int i = 1; i < nx - 1; i++)
            {
                onsite_G[i] = new DoubleComplexMatrix(ny,ny);
                hopping_G[i] = new DoubleComplexMatrix(ny,ny);
            }

            onsite_G[0] = Calculate_Boundary_Conditions(Boundary.left, energy);
            onsite_G[nx - 1] = Calculate_Boundary_Conditions(Boundary.right, energy);
            LHSboundary = onsite_G[0];
            RHSboundary = onsite_G[nx - 1];
        }

        void Add_Slice(double energy, int slice_no)
        {
            DoubleComplexMatrix new_slice = Create_New_Slice(energy, slice_no);
            Bind_Slice(slice_no, new_slice);
        }

        private void Bind_Slice(int slice_no, DoubleComplexMatrix slice_to_bind)
        {
            DoubleComplexMatrix new_hopping = Generate_Hopping_Hamiltonian();
            DoubleComplexMatrix new_hopping_trans = Generate_Hopping_Hamiltonian().Transpose();

            DoubleComplexMatrix temp_matrix = new DoubleComplexMatrix(slice_to_bind - Product(Product(new_hopping, onsite_G[slice_no - 1]), new_hopping_trans));

            // factorise the matrix
            DoubleComplexLUFact factorisation = new DoubleComplexLUFact(temp_matrix);

            // the inverse of this matrix is the onsite Green's function for the new slice
            onsite_G[slice_no] = factorisation.Inverse();

            // the first hop from n to n + 1 uses the onsite Green's function as G_(i,n+1) = G_(i,n) H' G(n+1,n+1) and when i = n, G_(i,n) = G(i,i)
            // must do this here before updating the onsite Green's functions
            hopping_G[slice_no - 1] = onsite_G[slice_no - 1];

            // update the old onsite Green's functions (not including the boundary slice)
            for (int i = 0; i < slice_no; i++)
                onsite_G[i] = onsite_G[i] + Product(Product(Product(Product(hopping_G[i], new_hopping), onsite_G[slice_no]), new_hopping_trans), hopping_G[i].Transpose());

            // calculate the new hopping Green's functions
            for (int i = 0; i < slice_no; i++)
                hopping_G[i] = Product(Product(hopping_G[i], new_hopping), onsite_G[slice_no]);
        }

        private DoubleComplexMatrix Create_New_Slice(double energy, int slice_no)
        {
            // work out where this slice is
            double x = (slice_no - (nx - 1) / 2) * dx;

            // generate the matrices for this slice
            DoubleComplexMatrix new_slice = Generate_Slice_Hamiltonian(x);

            return (energy * DoubleMatrix.Identity(ny)) - new_slice;
        }

        private DoubleComplexMatrix Calculate_Boundary_Conditions(Boundary boundary, double energy)
        {
            DoubleComplexMatrix bc_matrix = new DoubleComplexMatrix(2 * ny, 2 * ny);

            // generate Hamiltonians for the boundaries
            DoubleComplexMatrix tmp_slice;
            if (boundary == Boundary.left)
                tmp_slice = Generate_Slice_Hamiltonian(-0.5 * (nx - 1) * dx);
            else if (boundary == Boundary.right)
                tmp_slice = Generate_Slice_Hamiltonian(0.5 * (nx - 1) * dx);
            else throw new NotImplementedException();
            DoubleComplexMatrix tmp_hopping_trans = Generate_Hopping_Hamiltonian().Transpose();

            // create temporary matrix for top-left of BC matrix
            DoubleComplexMatrix tmp_matrix = Product(tmp_hopping_trans, (new DoubleComplexMatrix(energy * DoubleMatrix.Identity(ny)) - tmp_slice));

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
            coefficient_matrix = Product(Product(tmp_hopping_trans.Transpose(), eig_vecs), eig_vals);
            DoubleComplexLUFact lu_fact = new DoubleComplexLUFact(coefficient_matrix);
            return new DoubleComplexMatrix(Product(eig_vecs, lu_fact.Inverse()));

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
                if (eig_decomp.EigenValue(solution_no).Imag > 0 && boundary == Boundary.left)
                    return true;
                else if (eig_decomp.EigenValue(solution_no).Imag < 0 && boundary == Boundary.right)
                    return true;
                else if (eig_decomp.EigenValue(solution_no).Imag == 0)
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

        DoubleComplexMatrix Generate_Hopping_Hamiltonian()
        {
            DoubleComplexMatrix result = new DoubleComplexMatrix(ny,ny);
            for (int i = 0; i < ny; i++)
                result[i, i] = alphax;

            return result;
        }

        DoubleComplexMatrix Generate_Slice_Hamiltonian(double x)
        {
            DoubleComplexMatrix result = new DoubleComplexMatrix(ny,ny);
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

        double split_width = 800;
        double split_length = 800;
        double V_max = 5.0;
        double V_0 = 0.0;
        double Get_Potential(double x, double y)
        {
            //return V_0 / alphax_prime + V_max / alphax_prime * Convert.ToDouble((Math.Abs(x) < 0.5 * split_length-dx/2) && (Math.Abs(y) > 0.5 * split_width+dy/2));
            /* if (x > -150 && y > 0)
                 return 1000000000000;
             else return 0;*/
            return potential_xy[(int)((x+2500) / dx), (int)((y+1000) / dy)] / alphax_prime;
            //return 0;
            //return V_max / alphax_prime * Math.Exp(-1.0 * (x * x) / (2*1500*1500) );
        }

        DoubleComplexMatrix Product(DoubleHermitianMatrix A, DoubleHermitianMatrix B)
        {
            return new DoubleComplexMatrix(NMathFunctions.Product(MatrixFunctions.ToGeneralMatrix(A), MatrixFunctions.ToGeneralMatrix(B)));
        }

        DoubleComplexMatrix Product(DoubleComplexMatrix A, DoubleHermitianMatrix B)
        {
            return new DoubleComplexMatrix(NMathFunctions.Product(A, MatrixFunctions.ToGeneralMatrix(B)));
        }

        DoubleComplexMatrix Product(DoubleHermitianMatrix A, DoubleComplexMatrix B)
        {
            return new DoubleComplexMatrix(NMathFunctions.Product(MatrixFunctions.ToGeneralMatrix(A), B));
        }

        DoubleComplexMatrix Product(DoubleComplexMatrix A, DoubleComplexMatrix B)
        {
            return NMathFunctions.Product(A, B);
        }


        #region unused_methods
        public override void Get_ChargeDensity(Solver_Bases.Layers.ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override SpinResolved_Data Get_ChargeDensity(Solver_Bases.Layers.ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(Solver_Bases.Layers.ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override DoubleVector Get_EnergyLevels(Solver_Bases.Layers.ILayer[] layers, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }
        #endregion
    }
}
