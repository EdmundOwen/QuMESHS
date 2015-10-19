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
    class Strip : Slab
    {
        protected double dx, dy, half_width;
        protected double cos_theta, sin_theta;

        public Strip(double zmin, double zmax, double dx, double dy, double width, double theta)
            : base (zmin, zmax)
        {
            // convert theta into radians
            theta *= 2.0 * Math.PI / 360.0;

            this.dx = dx; this.dy = dy; this.half_width = 0.5 * width;
            this.cos_theta = Math.Cos(theta); this.sin_theta = Math.Sin(theta);

            // if the strip is parallel to one of the cardinal axes, set min/max values
            if (theta == 0.0 || theta == Math.PI)
            {
                xmin = -1.0 * half_width + dx; xmax = half_width + dx;
            }
            else if (theta == 0.5 * Math.PI || theta == 1.5 * Math.PI)
            {
                ymin = -1.0 * half_width + dy; ymax = half_width + dy;
            }
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = (x - dx) * cos_theta; double yprime = (y - dy) * sin_theta;
            // rotate the point by -theta and calculate if x is within half the width
            bool instrip = (Math.Abs(xprime + yprime) < half_width);
            return base.InLayer(0.0, 0.0, z) && instrip;
        }

        public new Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.strip; }
        }
    }
}
