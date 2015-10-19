/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

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
