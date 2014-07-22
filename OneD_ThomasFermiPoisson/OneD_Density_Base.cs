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
        public OneD_Density_Base(double temperature, double dz, int nz, double zmin)
            : base(temperature, 1.0, 1.0, dz, 1, 1, nz, 0.0, 0.0, zmin)
        { }

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

        public override double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z), temperature_input);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }

        public override void Close()
        {
            Console.WriteLine("Closing density solver");
        }
    }
}
