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
using CenterSpace.NMath.Core;
using Solver_Bases;
using Solver_Bases.Layers;

namespace TwoD_ThomasFermiPoisson
{
    public class TwoD_ThomasFermiSolver : TwoD_Density_Base
    {
        public TwoD_ThomasFermiSolver(Experiment exp)
            : this(exp, Carrier.electron)
        {
        }

        public TwoD_ThomasFermiSolver(IExperiment exp, Carrier carrier_type)
            : base(exp, carrier_type)
        {
        }

        public TwoD_ThomasFermiSolver(IExperiment exp, Plane plane, double plane_pos)
            : base(exp, plane, plane_pos)
        {
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // do not add anything to the density if on the edge of the domain
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    double x = dx * i + xmin;
                    double y = dy * j + ymin;
                    
                    // get the relevant layer
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, plane, x, y, pos_z);
                    
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
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    // leave the edges zeroed
                    if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
                        continue;

                    double x = dx * i + xmin;
                    double y = dy * j + ymin;

                    // get the relevant layer and if it's frozen out, don't recalculate the dopent charge
                    ILayer current_Layer = Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, plane, x, y, pos_z);

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
            Band_Data h_spin_summed = h.Spin_Summed_Data;

            System.IO.StreamWriter sw = new System.IO.StreamWriter(outfile);
            for (int i = 0; i < h_spin_summed.mat.Cols; i++)
                for (int j = 0; j < h_spin_summed.mat.Rows; j++)
                {
                    sw.Write(h_spin_summed.mat[j, i].ToString() + '\t');
                    if (j == h_spin_summed.mat.Rows - 1)
                        sw.WriteLine();
                }

            sw.Close();
        }

        public override void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot)
        {
            return new DoubleVector(1, double.NaN);
        }

        public Band_Data Get_KS_KE(ILayer[] layers, Band_Data chem_pot)
        {
            // return zeros (as Thomas-Fermi solution does not use wave functions)
            return new Band_Data(nx, ny, 0.0); ;
        }
    }
}
