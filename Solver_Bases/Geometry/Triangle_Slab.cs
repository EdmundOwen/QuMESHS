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
            this.x0 = x0; this.y0 = y0;
            this.theta1 = theta1; this.theta2 = theta2;
        }

        public override bool InLayer(double x, double y, double z)
        {
            // transform x and y so that the triangles point is at (0, 0)
            double xprime = x - x0; double yprime = y - y0;
            // and check that this point is within the triangle bounded from theta1 to theta2
            bool intriangle = (Math.Atan2(yprime, xprime) > theta1 && Math.Atan2(yprime, xprime) < theta2);
            return base.InLayer(z) && intriangle;
        }

        public virtual Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.triangle_slab; }
        }
    }
}
