using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Mixing_Schedulers
{
    class Constant_Scheduler : Scheduler_Base
    {
        double alpha;
        double zeta;

        public Constant_Scheduler(Dictionary<string, object> dict)
        {
            alpha = (double)dict["alpha"];
            zeta = (double)dict["zeta"];
        }

        public override double[] Get_Mixing_Parameter(int count)
        {
            return new double[] { alpha, zeta };
        }
    }
}
