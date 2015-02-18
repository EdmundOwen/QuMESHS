using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public abstract class FlexPDE_Base : Potential_Base
    {
        protected string output_file;

        public FlexPDE_Base(bool external_code, string external_input, string external_location, double tol)
            : base(external_code, external_input, external_location, tol) 
        {
            output_file = external_input;
            // for flexpde, we need the -S parameter included in the input
            base.external_input = "-S " + external_input;
        }

        protected override string[] Trim_Potential_File(string[] lines)
        {
            // work out where the data starts (this is flexPDE specific)
            int first_line = 0;
            for (int i = 0; i < lines.Length; i++)
                if (lines[i].StartsWith("}"))
                {
                    first_line = i + 1;
                    break;
                }

            // trim off the first lines which contain no data
            string[] data = new string[lines.Length - first_line];
            for (int i = first_line; i < lines.Length; i++)
                data[i - first_line] = lines[i];
            return data;
        }

        public abstract void Create_FlexPDE_File(double top_bc, double split_bc, double split_width, double surface, double bottom_bc, string output_file);
        public abstract void Create_NewtonStep_File(double top_bc, double split_bc, double split_width, double surface, double bottom_bc, string output_file, double t);
    }
}
