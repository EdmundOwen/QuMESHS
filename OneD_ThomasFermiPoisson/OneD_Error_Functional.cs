using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases;

namespace OneD_ThomasFermiPoisson
{
    class OneD_Error_Functional : DoubleFunctional
    {
        OneD_PoissonSolver pot_solv;
        OneD_ThomasFermiSolver dens_solv;

        int nz;
        int count = 0;

        public OneD_Error_Functional(OneD_ThomasFermiSolver dens_solver, OneD_PoissonSolver pot_solver, int nz)
            : base(nz)
        {
            this.dens_solv = dens_solver;
            this.pot_solv = pot_solver;
            this.nz = nz;
        }

        /// <summary>
        /// returns the error between the given density and a new density calculated from the potential due to this density
        /// </summary>
        public override double Evaluate(DoubleVector density)
        {
            // calculate the density for the potential given by the old density
            SpinResolved_DoubleVector new_density = dens_solv.Get_OneD_Density(pot_solv.Get_Potential(density));
            double error = (new_density.Spin_Summed_Vector - density).TwoNormSquared();

            Console.WriteLine("Iteration:\t" + count.ToString() + "\tError:\t" + error.ToString());
            count++;

            return error;
        }

        public override DoubleVector Gradient(DoubleVector x)
        {
            return base.Gradient(x);
        }

        public override double Laplacian(DoubleVector x)
        {
            throw new NotImplementedException();
        }
    }
}
