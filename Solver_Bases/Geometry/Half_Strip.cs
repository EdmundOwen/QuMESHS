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
    class Half_Strip : Strip
    {
        public Half_Strip(double zmin, double zmax, double dx, double dy, double width, double theta)
            : base (zmin, zmax, dx, dy, width, theta)
        {
            // convert theta into radians
            theta *= 2.0 * Math.PI / 360.0;

            // if the half-strip is parallel to one of the cardinal axes, set min/max values for the end section
            if (theta == 0.0)
                ymax = dy;
            else if (theta == Math.PI)
                ymin = dy;
            else if (theta == 0.5 * Math.PI)
                xmax = dx;
            else if (theta == 1.5 * Math.PI)
                xmin = dx;
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = (x - dx) * sin_theta; double yprime = (y - dy) * cos_theta;
            // rotate the point by -theta? and calculate if x is in the correct half-slab of the strip 
            // (the base class calculates whether it's in the strip or not)
            bool inhalfslab = (xprime + yprime < 0);
            return base.InLayer(x, y, z) && inhalfslab;
        }

        public new Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.half_strip; }
        }
    }
}
