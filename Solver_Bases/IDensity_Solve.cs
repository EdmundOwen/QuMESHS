using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;

namespace Solver_Bases
{
    public interface IDensity_Solve
    {
        SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot);
        void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data carrier_density, ref SpinResolved_Data dopent_density, Band_Data chem_pot);
        void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
        SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot);
        void Output(SpinResolved_Data data, string filename);

        CenterSpace.NMath.Core.DoubleVector Get_EnergyLevels(ILayer[] layers, Band_Data chem_pot);

        void Update_DFT_Potential(SpinResolved_Data car_dens);
        void Reset_DFT_Potential();
        void Print_DFT_diff(SpinResolved_Data car_dens);
        Band_Data DFT_diff(SpinResolved_Data car_dens);
        double DFT_Mixing_Parameter { get; set;  }
        Band_Data DFT_Potential { get; }

        double Get_XC_Potential(double charge_density);
        Band_Data Get_XC_Potential(SpinResolved_Data charge_density);
        Band_Data Get_XC_Potential_Deriv(SpinResolved_Data charge_density);

        double Unit_Charge { get; }
        double Mass { get; }

        void Close();
    }
}
