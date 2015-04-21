using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    class Metal_Layer : Layer
    {
        public Metal_Layer(IGeom geom, int layer_no)
            : base(geom, layer_no)
        { }

        protected override void Set_Material_Parameters()
        {
            material = Material.Metal;
            permitivity = 0.0;

            // no band gap in a metal...
            this.band_gap = 0.0;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
        }
    }
}
