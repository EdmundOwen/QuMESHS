using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Mixing_Schedulers
{
    abstract class Scheduler_Base : IScheduler
    {
        public double Get_Mixing_Parameter(int count, int type)
        {
            return Get_Mixing_Parameter(count)[type];
        }

        public double Get_Mixing_Parameter(int count, string type)
        {
            int type_index = -1;
            if (type == "alpha")
                type_index = 0;
            else if (type == "zeta")
                type_index = 1;
            else
                throw new NotImplementedException("Error - Mixing parameter of type \"" + type + "\" is not currently implemented");

            return Get_Mixing_Parameter(count)[type_index];
        }

        public abstract double[] Get_Mixing_Parameter(int count);
    }
}
