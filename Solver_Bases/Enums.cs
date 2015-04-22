using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public enum Spin
    {
        Up = 0,
        Down = 1
    }

    public enum Direction
    {
        x = 0,
        y = 1,
        z = 2
    }

    public enum Plane
    {
        xy = 0,
        yz = 1,
        zx = 2
    }

    public enum Material
    {
        GaAs,
        Al03GaAs,
        AlGaAs,
        AlAs,
        In075GaAs,
        InAlAs,
        InGaAs,
        InAs,
        PMMA,
        Air,
        Metal,
        Substrate
    }

    public enum Dopent
    {
        donor,
        acceptor
    }

    public enum Geometry_Type
    {
        triangle_slab,
        half_slab,
        strip,
        half_strip,
        slab,
        sheet, 
        composite
    }
}
