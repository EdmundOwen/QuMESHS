using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public class PMMA_Layer : Layer
    {
        public PMMA_Layer(IGeom geom, int layer_no)
            : base(geom, layer_no)
        { }

        protected override void Set_Material_Parameters()
        {
            material = Material.PMMA;
            permitivity = Physics_Base.epsilon_pmma * Physics_Base.epsilon_0;

            // set the PMMA band gap
            this.band_gap = 4400.0;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
        }
    }
}
