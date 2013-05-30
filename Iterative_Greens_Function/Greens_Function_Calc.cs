using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OpenCLNet;

namespace Iterative_Greens_Function
{
    class Greens_Function_Calc
    {
        int no_Energy_Count, no_Slices;
        double dE;
        Mem hop_mat, hop_mat_conj, H0_slice;

        int slice_width;

        public void Initialise()
        {
        }

        public Mem[] Iterate(out Mem G_0n, out Mem G_n0, Mem Potential, double Energy)
        {
            Mem[] G_ii = new Mem[no_Slices];
            Mem[] G_in = new Mem[no_Slices];
            Mem[] G_ni = new Mem[no_Slices];

            // Calculate Green's function boundary conditions

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

        void Iterate_Slice(ref Mem[] G_ii, ref Mem[] G_in, ref Mem[] G_ni, Mem Potential, int slice_no, double Energy)
        {
            // Calculate G_n+1,n+1 (MacKinnon 5a)
            Mem E_minus_H = Calculate_Hamiltonian(Potential, slice_no, Energy);
            G_ii[slice_no] = Calculate_New_Greens_Function(E_minus_H, G_ii[slice_no - 1]);
            
            // Calculate influence of the new slice on the previously calculated results
            Calculate_Connections(G_ii[slice_no], ref G_ii, ref G_in, ref G_ni, Potential, slice_no, Energy);
        }
        
        /// <summary>
        /// Calculates the Green's functions from each slice to the end and 
        /// how the new diagonal Green's functions are influenced by the new
        /// slice
        /// </summary>
        void Calculate_Connections(Mem G_nplusnplus, ref Mem[] G_ii, ref Mem[] G_in, ref Mem[] G_ni, Mem Potential, int slice_no, double Energy)
        {
            // Calculate G_n+1,i 
            for (int i = 0; i < slice_no; i++)
                G_ni[i] = ViennaCL.TripMatProd(G_ii[slice_no], hop_mat_conj, G_ni[i]);   // (MacKinnon 5c)

            // Calculate G_i,i ( note that G_ni has been updated to include the next slice so
            // it is technically G_{n+1, i} )
            for (int i = 1; i < slice_no; i++)
                G_ii[i] = ViennalCL.Add(G_ii[i], ViennalCL.TripProd(G_in[i], hop_mat, G_ni[i]));      // (MacKinnon 5d with i = j)
                    
            // Calculate G_i,n+1
            for (int i = 0; i < slice_no + 1; i++)
                G_in[i] = ViennaCL.TripProd(G_in[i], hop_mat, G_ii[slice_no]);  // (MacKinnon 5d)
        }

        Mem Calculate_Hamiltonian(Mem Potential, int slice_no, double Energy)
        {
            int slice_start = slice_no * slice_width;

            // Get potential from memory
            Mem pot_Slice = ViennaCL.GetSlice(Potential, slice_start, slice_start + slice_width, 1);
            pot_Slice = ViennaCL.Subtract(ViennaCL.Prod(Energy, ViennaCL.UnitVector), pot_Slice);

            // Construct (E - V) - H_0
            return ViennaCL.Subtract(ViennaCL.Prod(ViennaCL.Identity, pot_Slice), H0_slice);
        }

        Mem Calculate_New_Greens_Function(Mem E_minus_H, Mem G_nn)
        {
            // Calculate V G_nn V\dag
            Mem result = ViennaCL.TripProd(hop_mat, G_nn, hop_mat_conj);
            result = ViennaCL.Add(E_minus_H, result);

            // Invert matrix
            return ViennaCL.Invert(result);
        }

        public double Calculate_Currents()
        {
            throw new NotImplementedException();
        }
    }
}
