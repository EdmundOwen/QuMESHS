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
