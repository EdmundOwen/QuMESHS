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
            exp.Initialise(25.0, 0.1, 0.0000001, 40);
            exp.Run();
        }
    }
}
