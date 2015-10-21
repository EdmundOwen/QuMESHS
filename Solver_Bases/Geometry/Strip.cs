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
