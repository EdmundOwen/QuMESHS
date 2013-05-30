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

            // Add output lead
            Calculate_Connections(G_ii[no_Slices - 1], ref G_ii, ref G_in, ref G_ni, Potential, j, Energy);

            // copy G_0n and G_n0 to the signiture inputs for current calculations
            G_0n = G_in[0];
            G_n0 = G_ni[0];

            return G_ii;
        }

        void Iterate_Slice(ref Mem[] G_ii, ref Mem[] G_in, ref Mem[] G_ni, Mem Potential, int slice_no, double Energy)
        {
            // Calculate G_n+1,n+1
            
            // Calculate influence of the new slice on the previously calculated results
            Calculate_Connections(G_nplusnplus, ref G_ii, ref G_in, ref G_ni, Potential, slice_no, Energy);
        }

        /// <summary>
        /// Calculates the Green's functions from each slice to the end and 
        /// how the new diagonal Green's functions are influenced by the new
        /// slice
        /// </summary>
        void Calculate_Connections(Mem G_nplusnplus, ref Mem[] G_ii, ref Mem[] G_in, ref Mem[] G_ni, Mem Potential, int slice_no, double Energy)
        {
            // Calculate G_i,n+1

            // Calculate G_i,i

            // Calculate G_n+1,i
        }

        public double Calculate_Currents()
        {
            throw new NotImplementedException();
        }
    }
}
