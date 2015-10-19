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
    public interface IPoisson_Solve
    {
        void Initiate_Poisson_Solver(Dictionary<string, double> device_dimensions, Dictionary<string, double> boundary_conditions);
        Band_Data Get_Potential(Band_Data density);
        Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi);
        Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data carrier_density, Band_Data dft_pot, Band_Data dft_calc);

        Band_Data Chemical_Potential { get; }
        double T { get; set; }
    }
}
