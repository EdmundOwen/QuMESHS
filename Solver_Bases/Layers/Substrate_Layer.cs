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
            Material specific_material = Material.AlGaAs;

            if (specific_material == Material.GaAs)
            {
                permitivity = Physics_Base.epsilon_r_GaAs * Physics_Base.epsilon_0;

                // set the GaAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
                this.band_gap = 1420.0;
                allow_donors = true;
                this.acceptor_energy = -680.0; this.donor_energy = 704.0;
            }
            else if (specific_material == Material.AlGaAs)
            {
                permitivity = Physics_Base.epsilon_r_AlGaAs * Physics_Base.epsilon_0;

                // set the AlGaAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
                this.band_gap = 1800.0;
                allow_donors = true;
                this.acceptor_energy = -859.0; this.donor_energy = 869.0;
            }
            else
                throw new NotImplementedException();
        }

        internal override void Set_Freeze_Out_Temperature()
        {
            // and set a default freeze-out temperature for the dopents of 70K
            freeze_out_temp = 70.0;
        }
    }
}
