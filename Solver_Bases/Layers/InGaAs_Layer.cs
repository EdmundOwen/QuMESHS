using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public class InGaAs_Layer : Layer
    {
        double x;
        double permitivity_bowing_ratio = 13.5;

        public InGaAs_Layer(IGeom geom, int layer_no, double alloy_ratio)
            : base(geom, layer_no)
        {
            this.x = alloy_ratio;
            // re-set material parameters with correct alloy ratio
            Set_Material_Parameters();
        }

        protected override void Set_Material_Parameters()
        {
            material = Material.InGaAs;
            permitivity = (x * Physics_Base.epsilon_r_InAs + (1 - x) * Physics_Base.epsilon_r_GaAs + permitivity_bowing_ratio * x * (1 - x)) * Physics_Base.epsilon_0;

            // set the InGaAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
            this.band_gap = (1519.2 - 1583.7 * x + 475.0 * x * x);
            allow_donors = false;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
            // and set a default freeze-out temperature for the dopents of 70K
            freeze_out_temp = 70.0;
        }
    }
}
