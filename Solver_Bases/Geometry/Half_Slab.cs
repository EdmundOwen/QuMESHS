using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    class Half_Slab : Slab
    {
        double dy;
        double cos_theta, sin_theta;

        public Half_Slab(double zmin, double zmax, double dy, double theta)
            : base (zmin, zmax)
        {
            this.dy = dy;
            this.cos_theta = Math.Cos(theta); this.sin_theta = Math.Sin(theta);
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = x; double yprime = y - dy;
            // rotate the point by -theta and calculate if x is negative
            bool inhalfslab = (xprime * cos_theta - yprime * sin_theta < 0);
            return base.InLayer(z) && inhalfslab;
        }

        public virtual Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.half_slab; }
        }
    }
}
