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
    class Half_Slab : Slab
    {
        double dx, dy;
        double cos_theta, sin_theta;

        public Half_Slab(double zmin, double zmax, double dx, double dy, double theta)
            : base (zmin, zmax)
        {
            // convert theta to radians
            theta *= 2.0 * Math.PI / 360.0;

            this.dx = dx; this.dy = dy;
            this.cos_theta = Math.Cos(theta); this.sin_theta = Math.Sin(theta);

            // if the slab is parallel to one of the cardinal axes, set min/max values
            if (theta == 0.0)
                xmin = dx;
            else if (theta == 0.5 * Math.PI)
                ymax = dy;
            else if (theta == Math.PI)
                xmax = dx;
            else if (theta == 1.5 * Math.PI)
                ymin = dy;
        }

        public Half_Slab(double zmin, double zmax, double dx, double theta)
            : this (zmin, zmax, dx, 0.0, theta)
        {
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = (x - dx) * cos_theta; double yprime = (y - dy) * sin_theta;
            // rotate the point by -theta and calculate if x is negative
            bool inhalfslab = (xprime + yprime < 0);
            return base.InLayer(0.0, 0.0, z) && inhalfslab;
        }

        public override Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.half_slab; }
        }
    }
}
