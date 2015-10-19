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
