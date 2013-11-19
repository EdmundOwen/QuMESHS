using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Geometry
{
    public interface IGeom
    {
        bool InLayer(double z);
        bool InLayer(double y, double z);
        bool InLayer(double x, double y, double z);
        double Xmin { get; } double Xmax { get; }
        double Ymin { get; } double Ymax { get; }
        double Zmin { get; } double Zmax { get; }
        Geometry_Type Get_Geometry { get; }
    }
}
