using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;

namespace Solver_Bases
{
    public interface IExperiment
    {
        void Run();
        void Initialise(Dictionary<string, object> input_dict);

        SpinResolved_Data Carrier_Density { get; }
        SpinResolved_Data Dopent_Density { get; }
        Band_Data Chemical_Potential { get; }
        ILayer[] Layers { get; }

        double Temperature { get; }

        int Nx_Dens { get; }
        double Dx_Dens { get; }
        double Xmin_Dens { get; }
        int Ny_Dens { get; }
        double Dy_Dens { get; }
        double Ymin_Dens { get; }
        int Nz_Dens { get; }
        double Dz_Dens { get; }
        double Zmin_Dens { get; }

        int Nx_Pot { get; }
        double Dx_Pot { get; }
        double Xmin_Pot { get; }
        int Ny_Pot { get; }
        double Dy_Pot { get; }
        double Ymin_Pot { get; }
        int Nz_Pot { get; }
        double Dz_Pot { get; }
        double Zmin_Pot { get; }
    }
}
