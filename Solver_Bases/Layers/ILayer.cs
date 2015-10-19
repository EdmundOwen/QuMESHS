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
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public interface ILayer
    {
        bool Dopents_Frozen_Out(double temperature);
        void Set_Dopents(double acceptor_concentration, double donor_concentration);

        #region Access_Geometry
        bool InLayer(double z);
        bool InLayer(double y, double z);
        bool InLayer(double x, double y, double z);
        Geometry_Type Geometry { get; }
        ILayer Get_Component(int component_no);
        int No_Components { get; }

        double Xmin { get; }
        double Xmax { get; }
        double Ymin { get; }
        double Ymax { get; }
        double Zmin { get; }
        double Zmax { get; }
        #endregion

        #region Layer_Properties
        int Layer_No { get; }
        Material Material { get; }
        double Permitivity { get; }
        double Band_Gap { get; }
        double Donor_Energy { get; }
        double Acceptor_Energy { get; }
        double Donor_Conc { get; }
        double Acceptor_Conc { get; }
        double Dopent_FreezeOut_T { get; }
        #endregion

        double Electron_Mass { get; }
        double Hole_Mass { get; }
    }
}
