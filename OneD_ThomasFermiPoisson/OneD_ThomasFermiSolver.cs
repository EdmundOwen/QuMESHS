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
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Analysis;
using Solver_Bases;
using Solver_Bases.Layers;
using System.Web;

namespace OneD_ThomasFermiPoisson
{
    public class OneD_ThomasFermiSolver : OneD_Density_Base
    {
        public OneD_ThomasFermiSolver(IExperiment exp, double dz, double zmin, int nz)
            : base(exp)
        {
            this.dz = dz;
            this.zmin = zmin;
            this.nz = nz;
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_density, ref SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            Band_Data car_dens_spin_summed = carrier_density.Spin_Summed_Data;
            Band_Data dop_dens_spin_summed = dopent_density.Spin_Summed_Data;

            for (int i = 0; i < nz; i++)
            {
                double z = dz * i + zmin;

                double local_dopent_density;
                double local_carrier_density = car_dens_spin_summed.vec[i];

                // get the relevant layer and if it's frozen out, don't recalculate the charge
                ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z);

                // calculate the density at the given point
                ZeroD_Density charge_calc = new ZeroD_Density(current_Layer, temperature);
                 if (!current_Layer.Dopents_Frozen_Out(temperature))
                 {
                    local_dopent_density = charge_calc.Get_DopentDensity(chem_pot.vec[i]);
                    local_carrier_density = charge_calc.Get_CarrierDensity(chem_pot.vec[i]);
                }
                else
                {
                    // if the density is frozen out, on the first step, this will add the carrier density to the dopent density
                    // to give a total, frozen-out charge.  After that, the local carrier density is set to zero and so this value
                    // should not change
                    local_dopent_density = dop_dens_spin_summed.vec[i];
                    local_carrier_density = charge_calc.Get_CarrierDensity(chem_pot.vec[i]);
                }
 
                 // as there is no spin dependence in this problem yet, just divide the charge into spin-up and spin-down components equally
                dopent_density.Spin_Down.vec[i] = 0.5 * local_dopent_density;
                dopent_density.Spin_Up.vec[i] = 0.5 * local_dopent_density;
                carrier_density.Spin_Down.vec[i] = 0.5 * local_carrier_density;
                carrier_density.Spin_Up.vec[i] = 0.5 * local_carrier_density;
            }
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            SpinResolved_Data tmp_car = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));
            SpinResolved_Data tmp_dop = new SpinResolved_Data(new Band_Data(new DoubleVector(nz)), new Band_Data(new DoubleVector(nz)));
            Get_ChargeDensity(layers, ref tmp_car, ref tmp_dop, chem_pot);

            for (int i = 0; i < nz; i++)
            {
                density.Spin_Up.vec[i] = tmp_car.Spin_Up.vec[i] + tmp_dop.Spin_Up.vec[i];
                density.Spin_Down.vec[i] = tmp_car.Spin_Down.vec[i] + tmp_dop.Spin_Down.vec[i];
            }
        }

        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            Get_ChargeDensity(layers, ref carrier_density, ref dopent_density, chem_pot);
            return carrier_density + dopent_density;
        }

        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }
    }
}
