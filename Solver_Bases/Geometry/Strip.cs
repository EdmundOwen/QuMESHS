using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    class Strip : Slab
    {
        double dy, half_width;
        double cos_theta, sin_theta;

        public Strip(double zmin, double zmax, double dy, double width, double theta)
            : base (zmin, zmax)
        {
            this.dy = dy; this.half_width = 0.5 * width;
            this.cos_theta = Math.Cos(theta); this.sin_theta = Math.Sin(theta);
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = x; double yprime = y - dy;
            // rotate the point by -theta and calculate if x is within half the width
            bool inhalfslab = (Math.Abs(xprime * cos_theta - yprime * sin_theta) < 0);
            return base.InLayer(z) && inhalfslab;
        }

        public new Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.strip; }
        }
    }
}
