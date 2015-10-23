/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
