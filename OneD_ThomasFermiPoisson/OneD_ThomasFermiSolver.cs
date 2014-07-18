using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using CenterSpace.NMath.Analysis;
using Solver_Bases;
using Solver_Bases.Layers;

namespace OneD_ThomasFermiPoisson
{
    public class OneD_ThomasFermiSolver : OneD_Density_Base
    {
        public OneD_ThomasFermiSolver(double temperature, double dz, int nz, double zmin)
            : base(temperature, dz, nz, zmin)
        {
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_density, ref SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            for (int i = 0; i < nz; i++)
            {
                double z = dz * i + zmin;

                double local_dopent_density;
                double local_carrier_density = carrier_density.Spin_Summed_Data.vec[i];

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
                    local_dopent_density = dopent_density.Spin_Summed_Data.vec[i];
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
    }
}
