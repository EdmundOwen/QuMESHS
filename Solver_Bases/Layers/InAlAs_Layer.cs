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
    public class InAlAs_Layer : Layer
    {
        double x;
        double permitivity_bowing_ratio = 12.5;

        public InAlAs_Layer(IGeom geom, int layer_no, double alloy_ratio)
            : base(geom, layer_no)
        {
            // if alloy ratio is below 53% In, this is an indirect band gap material and I don't have the parameters for it
            if (alloy_ratio < 0.53)
                throw new NotImplementedException("Error - Below an indium concentration of 53%, the material has an indirect band gap.\nMaterial properties are not implemented for this material");

            this.x = alloy_ratio;
            // re-set material parameters with correct alloy ratio
            Set_Material_Parameters();
        }

        protected override void Set_Material_Parameters()
        {
            material = Material.InAlAs;
            permitivity = (x * Physics_Base.epsilon_r_InAs + (1 - x) * Physics_Base.epsilon_r_AlAs + permitivity_bowing_ratio * x * (1 - x)) * Physics_Base.epsilon_0;

            // set the InAlAs band gap and acceptor/donor energies are positivie and show how far from the band gap centre the donors are
            this.band_gap = 2640.0 - 2280.0 * x;
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
