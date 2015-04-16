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
            Material specific_material = Material.GaAs;
            ILayer tmp_layer;

            if (specific_material == Material.GaAs)
                tmp_layer = new GaAs_Layer(new Slab(0.0, 0.0), -1);
            else if (specific_material == Material.AlGaAs)
                // create a temporary gaas layer to obtain the material properties from (assume typical 33% AlGaAs)
                tmp_layer = new AlGaAs_Layer(new Slab(0.0, 0.0), -1, 0.33);
            else
                throw new NotImplementedException();

            permitivity = tmp_layer.Permitivity;

            // set the GaAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
            this.band_gap = tmp_layer.Band_Gap;
            allow_donors = true;
            this.acceptor_energy = tmp_layer.Acceptor_Energy; this.donor_energy = tmp_layer.Donor_Energy;
        }

        internal override void Set_Freeze_Out_Temperature()
        {
            // and set a default freeze-out temperature for the dopents of 70K
            freeze_out_temp = 70.0;
        }
    }
}
