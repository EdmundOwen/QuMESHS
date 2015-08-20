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
