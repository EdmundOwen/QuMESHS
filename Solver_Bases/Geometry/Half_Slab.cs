using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    class Half_Slab : Slab
    {
        double dx;
        double cos_theta, sin_theta;

        public Half_Slab(double zmin, double zmax, double dx, double dy, double theta)
            : this (zmin, zmax, dx - dy / Math.Tan(theta), theta)
        {
        }

        public Half_Slab(double zmin, double zmax, double dx, double theta)
            : base (zmin, zmax)
        {
            this.dx = dx;
            this.cos_theta = Math.Cos(theta); this.sin_theta = Math.Sin(theta);
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform y coordinate
            double xprime = x - dx; double yprime = y;
            // rotate the point by -theta and calculate if x is negative
            bool inhalfslab = (xprime * cos_theta - yprime * sin_theta < 0);
            return base.InLayer(z) && inhalfslab;
        }

        public new Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.half_slab; }
        }
    }
}
