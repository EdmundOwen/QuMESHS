using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Mixing_Schedulers
{
    public interface IScheduler
    {
        double[] Get_Mixing_Parameter(int count);
        double Get_Mixing_Parameter(int count, int type);
        double Get_Mixing_Parameter(int count, string type);
    }
}
