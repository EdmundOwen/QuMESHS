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
