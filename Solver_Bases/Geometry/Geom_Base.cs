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

namespace Solver_Bases.Geometry
{
    public abstract class Geom_Base
    {
        protected double xmin = -1.0 * double.MaxValue;
        protected double ymin = -1.0 * double.MaxValue;
        protected double zmin = -1.0 * double.MaxValue;
        protected double xmax = double.MaxValue;
        protected double ymax = double.MaxValue;
        protected double zmax = double.MaxValue;

        public double Xmin { get { return xmin; } } public double Xmax { get { return xmax; } }
        public double Ymin { get { return ymin; } } public double Ymax { get { return ymax; } }
        public double Zmin { get { return zmin; } } public double Zmax { get { return zmax; } }
    }
}
