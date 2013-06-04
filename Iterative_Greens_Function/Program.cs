using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Iterative_Greens_Function
{
    class Program
    {
        static void Main(string[] args)
        {
            Dictionary<string, object> input_dict = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref input_dict, "Input_Parameters.txt");

            Experiment exp = new Experiment();
            exp.Initialise_Experiment(input_dict);
            exp.Run();
        }
    }
}
