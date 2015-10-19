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
using Solver_Bases.Layers;

namespace Solver_Bases
{
    public interface IDensity_Solve
    {
        SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot);
        void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_density, ref SpinResolved_Data dopent_density, Band_Data chem_pot);
        void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
        SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot);
        void Output(SpinResolved_Data data, string filename);

        CenterSpace.NMath.Core.DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot);

        void Update_DFT_Potential(SpinResolved_Data car_dens);
        void Reset_DFT_Potential();
        void Print_DFT_diff(SpinResolved_Data car_dens);
        Band_Data DFT_diff(SpinResolved_Data car_dens);
        double DFT_Mixing_Parameter { get; set;  }
        Band_Data DFT_Potential { get; }

        double Get_XC_Potential(double charge_density);
        Band_Data Get_XC_Potential(SpinResolved_Data charge_density);
        Band_Data Get_XC_Potential_Deriv(SpinResolved_Data charge_density);

        double Unit_Charge { get; }
        double Mass { get; }

        void Close();
    }
}
