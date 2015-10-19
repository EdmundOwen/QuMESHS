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
using System.Threading.Tasks;
using Solver_Bases;

namespace Iterative_Greens_Function_Test
{
    public class Experiment : Experiment_Base
    {
        public Experiment()
            : base()
        {
        }

        public override void Initialise(Dictionary<string, object> input_dict)
        {
            //base.Initialise(input_dict);
        }

        public override bool Run()
        {
            Iterative_Greens_Function iter = new Iterative_Greens_Function(this);
            iter.Iterate();

            //throw new NotImplementedException();
            return true;
        }
  

  /* public double[,] CalcDoS(double energy, double[,] potential)
        {
            Iterative_Greens_Function iter = new Iterative_Greens_Function(this, potential);

            double[,] DoS = iter.GetDoS(energy);

            return DoS;
        }*/

        protected override void Initialise_DataClasses(Dictionary<string, object> input_dict)
        {
            throw new NotImplementedException();
        }

        protected override void Initialise_from_Hot_Start(Dictionary<string, object> input_dict)
        {
            throw new NotImplementedException();
        }

        protected override bool Run_Iteration_Routine(IDensity_Solve dens_solv, IPoisson_Solve pois_solv, double tol, int max_iterations)
        {
            throw new NotImplementedException();
        }
    }
}
