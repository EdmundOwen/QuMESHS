using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    /// <summary>
    /// A horizontal slab
    /// </summary>
    class Slab : Geom_Base, IGeom
    {
        public Slab(double zmin, double zmax)
        {
            this.zmin = zmin; this.zmax = zmax;
        }

        public bool InLayer(double z)
        {
            return InLayer(0.0, 0.0, z);
        }

        public bool InLayer(double y, double z)
        {
            return InLayer(0.0, y, z);
        }

        public virtual bool InLayer(double x, double y, double z)
        {
            return (z <= zmax && z > zmin);
        }

        public Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.slab; }
        }
    }
}