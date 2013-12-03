using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public class Air_Layer : Layer
    {
        public Air_Layer(IGeom geom, int layer_no)
            : base(geom, layer_no)
        { }

        protected override void Set_Material_Parameters()
        {
            material = Material.Air;
            permitivity = Physics_Base.epsilon_0;

            // ok, so this is the actual band gap for air... may as well be infinite but that's the way I'm doing this
            this.band_gap = 500000000.0;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
        }
    }
}
