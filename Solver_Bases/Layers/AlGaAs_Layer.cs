using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public class AlGaAs_Layer : Layer
    {
        public AlGaAs_Layer(IGeom geom, int layer_no)
            : base(geom, layer_no)
        { }

        protected override void Set_Material_Parameters()
        {
            material = Material.AlGaAs;
            permitivity = Physics_Base.epsilon_r_AlGaAs * Physics_Base.epsilon_0;

            // set the AlGaAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
            this.band_gap = 1835.0;
            allow_donors = true;
            this.acceptor_energy = -877.0; this.donor_energy = 887.0;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
            // and set a default freeze-out temperature for the dopents of 70K
            freeze_out_temp = 70.0;
        }
    }
}
