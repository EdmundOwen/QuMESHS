﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public class AlGaAs_Layer : Layer
    {
        double x;
        double permitivity_bowing_ratio = 0.0;

        public AlGaAs_Layer(IGeom geom, int layer_no, double alloy_ratio)
            : base(geom, layer_no)
        {
            this.x = alloy_ratio;
            // re-set material parameters with correct alloy ratio
            Set_Material_Parameters();
        }

        protected override void Set_Material_Parameters()
        {
            material = Material.AlGaAs;
            permitivity = ((1 - x) * Physics_Base.epsilon_r_GaAs + x * Physics_Base.epsilon_r_AlAs + permitivity_bowing_ratio * x * (1 - x)) * Physics_Base.epsilon_0;

            // set the AlGaAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
            if (x < 0.45)
                this.band_gap = 1424.0 + 1247.0 * x;
            else
                this.band_gap = 1900.0 + 125.0 * x + 143.0 * x * x;

            allow_donors = true;
            this.acceptor_energy = -0.5 * band_gap + 35.0; this.donor_energy = 0.5 * band_gap - 5.0;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
            // and set a default freeze-out temperature for the dopents of 70K
            freeze_out_temp = 70.0;
        }
    }
}
