using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public abstract class FlexPDE_Base : Potential_Base
    {
        protected string flexpde_options;
        protected string flexpde_script;
        protected string flexpde_output;

        public FlexPDE_Base(bool external_code, Dictionary<string, object> input)
            : base(external_code) 
        {
            if (!input.ContainsKey("FlexPDE_location") || !input.ContainsKey("FlexPDE_file"))
                throw new Exception("Error - To use FlexPDE you must provide the location of FlexPDE in \"FlexPDE_location\" and the script location in \"FlexPDE_file\"");

            string flexpde_location = (string)input["FlexPDE_location"];

            this.initcalc_location = flexpde_location;
            this.newton_location = flexpde_location;

            this.flexpde_script = (string)input["FlexPDE_file"];
            this.flexpde_output = flexpde_script.Remove(flexpde_script.Length - 4) + (string)input["output_suffix"] + ".pg6";

            // for flexpde, we need the -S parameter included in the input
            flexpde_options = "-S " + flexpde_script;
        }

        protected override string[] Trim_Potential_File(string[] lines)
        {
            // copy the pg6 file to a saved location
            System.IO.File.Copy(flexpde_script.Remove(flexpde_script.Length - 4) + ".pg6", flexpde_output, true);

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

        protected override Band_Data Get_Pot_On_Regular_Grid(Band_Data density)
        {
            throw new NotImplementedException();
        }
    }
}
