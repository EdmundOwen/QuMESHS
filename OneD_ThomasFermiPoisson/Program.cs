using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace OneD_ThomasFermiPoisson
{
    class Program
    {
        static void Main(string[] args)
        {
            Experiment exp = new Experiment();
            exp.Initialise(25.0, 0.01, 0.00000000001, 10);
            exp.Run();
        }
    }
}
