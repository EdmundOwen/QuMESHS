using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    class Strip : Slab
    {
        double dx, dy, half_width;
        double cos_theta, sin_theta;

        public Strip(double zmin, double zmax, double dx, double dy, double width, double theta)
            : base (zmin, zmax)
        {
            // convert theta into radians
            theta *= 2.0 * Math.PI / 360.0;

            this.dx = dx; this.dy = dy; this.half_width = 0.5 * width;
            this.cos_theta = Math.Cos(theta); this.sin_theta = Math.Sin(theta);
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = (x - dx) * cos_theta; double yprime = (y - dy) * sin_theta;
            // rotate the point by -theta and calculate if x is within half the width
            bool inhalfslab = (Math.Abs(xprime + yprime) < 0);
            return base.InLayer(0.0, 0.0, z) && inhalfslab;
        }

        public new Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.strip; }
        }
    }
}
