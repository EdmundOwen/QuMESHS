using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    class Substrate_Layer : Layer
    {
        public Substrate_Layer(IGeom geom, int layer_no)
            : base(geom, layer_no)
        { }

        protected override void Set_Material_Parameters()
        {
            material = Material.Substrate;

            // set substrate parameters to undoped GaAs
            this.band_gap = 1420.0;
        }

        internal override void Set_Freeze_Out_Temperature() 
        { 
        }
    }
}
