/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public static class Physics_Base
    {
        // physical constants
        public const double hbar = 0.658211814;                 // (meV) (ps)
        public const double m_e = 5.68562958e-3;                // (meV) (ps)^2 (nm)^-2 of the free electron mass
        public const double q_e = 160.217646;                   // (zC) is positive as in these definitions it is the elementary charge
        public const double energy_V_to_meVpzC = 6.2415093;     // conversion factor from V to meV per zC
        public const double epsilon_0 = 1.41859713;             // (zC)^2 (nm)^-1 (meV)^-1 for vacuum
        public const double epsilon_r = 13.0;                   // relative permittivity for vacuum -> GaAs
        public const double epsilon_r_GaAs = 12.9;              // relative permittivity for vacuum -> GaAs
        public const double epsilon_r_InAs = 14.6;              // relative permittivity for vacuum -> InAs
        public const double epsilon_r_AlAs = 10.1;              // relative permittivity for vacuum -> AlAs
        public const double epsilon_pmma = 2.6;                 // relative permittivity for vacuum -> PMMA
        public const double epsilon_al2o3 = 9.3;                 // relative permittivity for vacuum -> Al2O3 (assuming epsilon_11)
        public const double epsilon = epsilon_0 * epsilon_r;    // (zC)^2 (nm)^-1 (meV)^-1 for GaAs
        public const double kB = 0.086173324;                   // (meV) (K)^-1

        /// <summary>
        /// gets 3D spin-resolved density of states for a given potential and energy
        /// </summary>
        public static double Get_3D_DensityofStates(double energy, double band_edge, Carrier carrier_type, double mass)
        {
            // if the energy is below the potential (i.e. below the conduction band) return zero
            if (carrier_type == Carrier.electron && energy < band_edge)
                return 0.0;
            else if (carrier_type == Carrier.hole && energy > band_edge)
                return 0.0;

            // generate density of states
            if (carrier_type == Carrier.electron)
                return Math.Pow(2.0 * mass * mass * mass * (energy - band_edge), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);
            else if (carrier_type == Carrier.hole)
                return Math.Pow(2.0 * mass * mass * mass * (band_edge - energy), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);
            else
                throw new NotImplementedException();
        }

    }
}
