using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Solver_Bases;

namespace Iterative_Greens_Function_Test
{
    class Experiment : Experiment_Base
    {
        public Experiment()
            : base()
        {
        }

        public override void Initialise(Dictionary<string, object> input_dict)
        {
            //base.Initialise(input_dict);
        }

        public override bool Run()
        {
            Iterative_Greens_Function iter = new Iterative_Greens_Function(this);
            iter.Iterate();

            throw new NotImplementedException();
        }


        protected override void Initialise_DataClasses(Dictionary<string, object> input_dict)
        {
            throw new NotImplementedException();
        }

        protected override void Initialise_from_Hot_Start(Dictionary<string, object> input_dict)
        {
            throw new NotImplementedException();
        }

        protected override bool Run_Iteration_Routine(IDensity_Solve dens_solv, IPoisson_Solve pois_solv, double tol, int max_iterations)
        {
            throw new NotImplementedException();
        }
    }
}
