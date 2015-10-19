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
    public interface IExperiment
    {
        bool Run();
        void Initialise(Dictionary<string, object> input_dict);

        SpinResolved_Data Carrier_Density { get; }
        SpinResolved_Data Dopent_Density { get; }
        Band_Data Chemical_Potential { get; }
        Band_Data GPhi { get; }
        Band_Data X { get; }
        ILayer[] Layers { get; }

        double Temperature { get; }
        double Current_Temperature { get; }

        int Nx_Dens { get; }
        double Dx_Dens { get; }
        double Xmin_Dens { get; }
        int Ny_Dens { get; }
        double Dy_Dens { get; }
        double Ymin_Dens { get; }
        int Nz_Dens { get; }
        double Dz_Dens { get; }
        double Zmin_Dens { get; }

        int Nx_Pot { get; }
        double Dx_Pot { get; }
        double Xmin_Pot { get; }
        int Ny_Pot { get; }
        double Dy_Pot { get; }
        double Ymin_Pot { get; }
        int Nz_Pot { get; }
        double Dz_Pot { get; }
        double Zmin_Pot { get; }
    }
}
