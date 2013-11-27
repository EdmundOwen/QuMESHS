using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;

namespace Solver_Bases
{
    public interface IDensity_Solve
    {
        void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
        void Output(SpinResolved_Data data, string filename);

        double Get_Chemical_Potential(double z, ILayer[] layers);
        double Get_Chemical_Potential(double z, ILayer[] layers, double temperature_input);
        double Get_Chemical_Potential(double y, double z, ILayer[] layers);
        double Get_Chemical_Potential(double y, double z, ILayer[] layers, double temperature_input);
        double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers);
        double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input);

        void Close();
    }
}
