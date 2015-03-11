using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public interface IPoisson_Solve
    {
        void Initiate_Poisson_Solver(Dictionary<string, double> device_dimensions, Dictionary<string, double> boundary_conditions);
        Band_Data Get_Chemical_Potential(Band_Data density);
        Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi);
        Band_Data Calculate_Newton_Step(SpinResolved_Data rho_prime, Band_Data gphi, SpinResolved_Data carrier_density, Band_Data dft_difference);
        Band_Data Calculate_Laplacian(Band_Data input_vec);

        Band_Data Chemical_Potential { get; }
        double T { get; set; }
    }
}
