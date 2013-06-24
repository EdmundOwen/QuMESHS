using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    class Program
    {
        static void Main(string[] args)
        {
            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");

            Experiment exp = new Experiment();
            exp.Initialise(inputs);
            exp.Run();


        }
    }
}
