using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases;

namespace TwoD_ThomasFermiPoisson
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Program starting");

            Console.WriteLine("Loading input parameters from file");
            Dictionary<string, object> inputs = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters.txt");
            Console.WriteLine("Input parameters loaded");

            Console.WriteLine("Performing density dopent calculation");
            Dictionary<string, object> inputs_init = new Dictionary<string, object>();
            Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs_init, "Input_Parameters_1D.txt");
            OneD_ThomasFermiPoisson.Experiment exp_init = new OneD_ThomasFermiPoisson.Experiment();
            exp_init.Initialise(inputs_init);
            exp_init.Run();
            Console.WriteLine("Calculated 1D density for dopents");

            Console.WriteLine("Starting experiment");
            Experiment exp = new Experiment();
            Console.WriteLine("Created experiment");
            exp.Initialise_Experiment(inputs);
            exp.Run();
            Console.WriteLine("Experiment complete");
        }
    }
}
