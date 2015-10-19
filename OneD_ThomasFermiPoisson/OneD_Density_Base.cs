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
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;

namespace OneD_ThomasFermiPoisson
{
    public abstract class OneD_Density_Base : Density_Base
    {
        IExperiment exp;
        protected double dz;
        protected double zmin;
        protected int nz;

        public OneD_Density_Base(IExperiment exp)
            : base(exp.Temperature)
        {
            this.exp = exp;

            this.dz = exp.Dz_Dens;
            this.zmin = exp.Zmin_Dens;
            this.nz = exp.Nz_Dens;
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            for (int i = 0; i < nz; i++)
            {
                double z = dz * i + zmin;

                // get the relevant layer and if it's frozen out, don't recalculate the dopent charge
                ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z);

                ZeroD_Density charge_calc = new ZeroD_Density(current_Layer, temperature);
                if (!current_Layer.Dopents_Frozen_Out(temperature))
                {
                    double local_dopent_density_deriv = charge_calc.Get_DopentDensityDeriv(chem_pot.vec[i]);
                    dopent_density_deriv.Spin_Up.vec[i] = 0.5 * local_dopent_density_deriv;
                    dopent_density_deriv.Spin_Down.vec[i] = 0.5 * local_dopent_density_deriv;
                }
                else
                {
                    dopent_density_deriv.Spin_Up.vec[i] = 0.0;
                    dopent_density_deriv.Spin_Down.vec[i] = 0.0;
                }

                carrier_density_deriv.Spin_Up.vec[i] = 0.5 * charge_calc.Get_CarrierDensityDeriv(chem_pot.vec[i]);
                carrier_density_deriv.Spin_Down.vec[i] = 0.5 * charge_calc.Get_CarrierDensityDeriv(chem_pot.vec[i]);
            }

            return carrier_density_deriv + dopent_density_deriv;
        }

        public double Get_Chemical_Potential(double z, ILayer[] layers, double temperature_input)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z), temperature_input);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }
    }
}
