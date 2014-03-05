using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;

namespace OneD_ThomasFermiPoisson
{
    class Program
    {
        static void Main(string[] args)
        {
            // set nmath license key

            Dictionary<string,object> input = new Dictionary<string,object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref input, "Input_Parameters.txt");

            Experiment exp = new Experiment();
            exp.Initialise(input);
            exp.Run();
        }
    }
}
