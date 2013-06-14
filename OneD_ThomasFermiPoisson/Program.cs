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
            exp.Initialise(5.0, 0.0000001, 0.00000000001, 50);
            exp.Run();
        }
    }
}
