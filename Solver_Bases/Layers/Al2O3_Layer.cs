using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    class Al2O3_Layer : Layer
    {
        public Al2O3_Layer(IGeom geom, int layer_no)
            : base(geom, layer_no)
        { }

        protected override void Set_Material_Parameters()
        {
            material = Material.Al2O3;
            permitivity = Physics_Base.epsilon_al2o3 * Physics_Base.epsilon_0;

            // set the PMMA band gap
            this.band_gap = 8900.0;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
        }
    }
}
