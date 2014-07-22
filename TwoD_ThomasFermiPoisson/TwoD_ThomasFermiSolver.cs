using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;
using Solver_Bases.Layers;

namespace TwoD_ThomasFermiPoisson
{
    public class TwoD_ThomasFermiSolver : TwoD_Density_Base
    {
        public TwoD_ThomasFermiSolver(Experiment exp)
            : base(exp)
        {
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                        continue;

                    double y = dy * i + ymin;
                    double z = dz * j + zmin;
                    
                    // get the relevant layer
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, y, z);
                    
                    // calculate the density at the given point
                    ZeroD_Density charge_calc = new ZeroD_Density(current_Layer, temperature);
                    double local_charge_density = charge_calc.Get_CarrierDensity(chem_pot.mat[i, j]);

                    // as there is no spin dependence in this problem yet, just divide the charge into spin-up and spin-down components equally
                    density.Spin_Down.mat[i, j] = 0.5 * local_charge_density;
                    density.Spin_Up.mat[i, j] = 0.5 * local_charge_density;
                }
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    // leave the edges zeroed
                    if (i == 0 || i == ny - 1 || j == 0 || j == nz - 1)
                        continue;

                    double y = dy * i + ymin;
                    double z = dz * j + zmin;

                    // get the relevant layer and if it's frozen out, don't recalculate the dopent charge
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, y, z);

                    ZeroD_Density charge_calc = new ZeroD_Density(current_Layer, temperature);
                    if (!current_Layer.Dopents_Frozen_Out(temperature))
                    {
                        double local_dopent_density_deriv = charge_calc.Get_DopentDensityDeriv(chem_pot.mat[i, j]);
                        dopent_density_deriv.Spin_Up.mat[i, j] = 0.5 * local_dopent_density_deriv;
                        dopent_density_deriv.Spin_Down.mat[i, j] = 0.5 * local_dopent_density_deriv;
                    }
                    else
                    {
                        dopent_density_deriv.Spin_Up.mat[i, j] = 0.0;
                        dopent_density_deriv.Spin_Down.mat[i, j] = 0.0;
                    }

                    carrier_density_deriv.Spin_Up.mat[i, j] = 0.5 * charge_calc.Get_CarrierDensityDeriv(chem_pot.mat[i, j]);
                    carrier_density_deriv.Spin_Down.mat[i, j] = 0.5 * charge_calc.Get_CarrierDensityDeriv(chem_pot.mat[i, j]);
                }

            return carrier_density_deriv + dopent_density_deriv;
        }

        public void Write_Out_Density(SpinResolved_Data h, string outfile)
        {
            System.IO.StreamWriter sw = new System.IO.StreamWriter(outfile);
            for (int i = 0; i < h.Spin_Summed_Data.mat.Cols; i++)
                for (int j = 0; j < h.Spin_Summed_Data.mat.Rows; j++)
                {
                    sw.Write(h.Spin_Summed_Data.mat[j, i].ToString() + '\t');
                    if (j == h.Spin_Summed_Data.mat.Rows - 1)
                        sw.WriteLine();
                }

            sw.Close();
        }

        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }
    }
}
