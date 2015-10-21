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
