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

        void Set_DFT_Potential(SpinResolved_Data car_dens);
        void Reset_DFT_Potential();
        void Print_DFT_diff(SpinResolved_Data car_dens);
        Band_Data DFT_diff(SpinResolved_Data car_dens);
        double DFT_Mixing_Parameter { get; set;  }

        double Get_Chemical_Potential(double z, ILayer[] layers);
        double Get_Chemical_Potential(double z, ILayer[] layers, double temperature_input);
        double Get_Chemical_Potential(double y, double z, ILayer[] layers);
        double Get_Chemical_Potential(double y, double z, ILayer[] layers, double temperature_input);
        double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers);
        double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input);

        void Close();
    }
}
