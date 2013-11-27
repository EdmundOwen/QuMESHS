using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    /// <summary>
    /// A special class of "Slab" which has a minimal width
    /// </summary>
    class Sheet : Geom_Base, IGeom
    {
        public Sheet(double z)
        {
            this.zmin = z - double.MinValue; this.zmax = z;
        }

        public bool InLayer(double z)
        {
            return InLayer(0.0, 0.0, z);
        }

        public bool InLayer(double y, double z)
        {
            return InLayer(0.0, y, z);
        }

        public bool InLayer(double x, double y, double z)
        {
            return (z <= zmax && z > zmin);
        }

        public Geometry_Type Get_Geometry
        {
            get { return Geometry_Type.sheet; }
        }
    }
}
