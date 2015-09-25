using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;

namespace OneD_ThomasFermiPoisson
{
    public class OneD_eh_DFTSolver : IDensity_Solve
    {
        OneD_DFTSolver electron_dens_calc;
        OneD_DFTSolver hole_dens_calc;

        public OneD_eh_DFTSolver(IExperiment exp)
        {
            electron_dens_calc = new OneD_DFTSolver(exp, Carrier.Electron);
            hole_dens_calc = new OneD_DFTSolver(exp, Carrier.Hole);
        }

        public void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data charge_density, Band_Data chem_pot)
        {
            SpinResolved_Data electron_density = charge_density.DeepenThisCopy();
            SpinResolved_Data hole_density = charge_density.DeepenThisCopy();

            // remove the electron density from the hole density SpinResolved_Data and vice versa
            for (int i = 0; i < electron_density.Spin_Summed_Data.Length; i++ )
            {
                if (electron_density.Spin_Summed_Data.vec[i] > 0.0)
                {
                    electron_density.Spin_Up[i] = 0.0;
                    electron_density.Spin_Down[i] = 0.0;
                }
                if (hole_density.Spin_Summed_Data.vec[i] < 0.0)
                {
                    hole_density.Spin_Up[i] = 0.0;
                    hole_density.Spin_Down[i] = 0.0;
                }
            }

            electron_dens_calc.Get_ChargeDensity(layers, ref electron_density, chem_pot);
            hole_dens_calc.Get_ChargeDensity(layers, ref hole_density, chem_pot);

            charge_density = electron_density + hole_density;
        }
         
        public SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_charge_density, SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            SpinResolved_Data electron_density = carrier_charge_density.DeepenThisCopy();
            SpinResolved_Data hole_density = carrier_charge_density.DeepenThisCopy();

            // remove the electron density from the hole density SpinResolved_Data and vice versa
            for (int i = 0; i < electron_density.Spin_Summed_Data.Length; i++)
            {
                if (electron_density.Spin_Summed_Data.vec[i] > 0.0)
                {
                    electron_density.Spin_Up[i] = 0.0;
                    electron_density.Spin_Down[i] = 0.0;
                }
                if (hole_density.Spin_Summed_Data.vec[i] < 0.0)
                {
                    hole_density.Spin_Up[i] = 0.0;
                    hole_density.Spin_Down[i] = 0.0;
                }
            }

            // must remove one lot of dopents due to double counting
            return electron_dens_calc.Get_ChargeDensity(layers, electron_density, dopent_charge_density, chem_pot) + hole_dens_calc.Get_ChargeDensity(layers, hole_density, dopent_charge_density, chem_pot) 
                    - dopent_charge_density;
        }

        public SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            SpinResolved_Data electron_density_deriv = carrier_density_deriv.DeepenThisCopy();
            SpinResolved_Data hole_density_deriv = carrier_density_deriv.DeepenThisCopy();

            // remove the electron density from the hole density SpinResolved_Data and vice versa
            for (int i = 0; i < electron_density_deriv.Spin_Summed_Data.Length; i++)
            {
                if (electron_density_deriv.Spin_Summed_Data.vec[i] > 0.0)
                {
                    electron_density_deriv.Spin_Up[i] = 0.0;
                    electron_density_deriv.Spin_Down[i] = 0.0;
                }
                if (hole_density_deriv.Spin_Summed_Data.vec[i] < 0.0)
                {
                    hole_density_deriv.Spin_Up[i] = 0.0;
                    hole_density_deriv.Spin_Down[i] = 0.0;
                }
            }

            return electron_dens_calc.Get_ChargeDensity_Deriv(layers, electron_density_deriv, dopent_density_deriv, chem_pot) + hole_dens_calc.Get_ChargeDensity_Deriv(layers, hole_density_deriv, dopent_density_deriv, chem_pot);
        }

        public double Zmin_Pot
        {
            set { electron_dens_calc.Zmin_Pot = value; hole_dens_calc.Zmin_Pot = value; }
        }
        public double Dz_Pot
        {
            set { electron_dens_calc.Dz_Pot = value; hole_dens_calc.Dz_Pot = value; }
        }

        /// <summary>
        /// for this method, it is assumed that the dopent density is frozen out (and therefore not altered) unless overriden
        /// </summary>
        public void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_charge_density, ref SpinResolved_Data dopent_charge_density, Band_Data chem_pot)
        {
            Get_ChargeDensity(layers, ref carrier_charge_density, chem_pot);
        }

        public void Output(SpinResolved_Data data, string filename)
        {
            throw new NotImplementedException();
        }

        public DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public void Set_DFT_Potential(SpinResolved_Data car_dens)
        {
        }

        public void Reset_DFT_Potential()
        {
        }

        public void Print_DFT_diff(SpinResolved_Data car_dens)
        {
        }

        public Band_Data DFT_diff(SpinResolved_Data car_dens)
        {
            throw new NotImplementedException();
        }

        public double DFT_Mixing_Parameter
        {
            get
            {
                return double.NaN;
            }
            set
            {
            }
        }

        public void Close()
        {
            electron_dens_calc.Close();
            hole_dens_calc.Close();
        }
    }
}
