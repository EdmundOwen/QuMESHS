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
    class Triangle_Slab : Slab
    {
        double x0, y0;
        double theta1, theta2;

        public Triangle_Slab(double zmin, double zmax, double x0, double y0, double theta1, double theta2)
            : base (zmin, zmax)
        {
            // convert thetas to radians
            theta1 *= 2.0 * Math.PI / 360.0;
            theta2 *= 2.0 * Math.PI / 360.0;

            this.x0 = x0; this.y0 = y0;
            this.theta1 = theta1; this.theta2 = theta2;
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform x and y so that the triangles point is at (0, 0)
            double xprime = x - x0; double yprime = y - y0;
            // and check that this point is within the triangle bounded from theta1 to theta2
            bool intriangle = (Math.Atan2(yprime, xprime) > theta1 && Math.Atan2(yprime, xprime) < theta2);
            return base.InLayer(0.0, 0.0, z) && intriangle;
        }

        public new Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.triangle_slab; }
        }
    }
}
