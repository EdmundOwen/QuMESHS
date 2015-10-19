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

using CenterSpace.NMath.Core;
using CenterSpace.NMath.Matrix;

namespace Iterative_Greens_Function_Test
{
    class Program
    {
        static void Main(string[] args)
        {
            // set nmath license key

            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "inputs.txt");

            Experiment exp = new Experiment();
            exp.Initialise(inputs);
            exp.Run();
        }
    }
}
