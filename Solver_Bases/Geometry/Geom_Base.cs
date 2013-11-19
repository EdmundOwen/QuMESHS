using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    public abstract class Geom_Base
    {
        protected double xmin = -1.0 * double.MaxValue;
        protected double ymin = -1.0 * double.MaxValue;
        protected double zmin = -1.0 * double.MaxValue;
        protected double xmax = double.MaxValue;
        protected double ymax = double.MaxValue;
        protected double zmax = double.MaxValue;

        public double Xmin { get { return xmin; } } public double Xmax { get { return xmax; } }
        public double Ymin { get { return ymin; } } public double Ymax { get { return ymax; } }
        public double Zmin { get { return zmin; } } public double Zmax { get { return zmax; } }
    }
}
