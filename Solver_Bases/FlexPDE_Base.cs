/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

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

        protected double pot_tol = 1e-6;
        protected double newton_tol = 1e-5;

        public FlexPDE_Base(bool external_code, Dictionary<string, object> input)
            : base(external_code) 
        {
            if (!input.ContainsKey("FlexPDE_location") || !input.ContainsKey("FlexPDE_file"))
                throw new Exception("Error - To use FlexPDE you must provide the location of FlexPDE in \"FlexPDE_location\" and the script location in \"FlexPDE_file\"");

            string flexpde_location = (string)input["FlexPDE_location"];

            this.initcalc_location = flexpde_location;
            this.newton_location = flexpde_location;

            // if an output suffix was not defined before, allocate it as a null string
            if (!input.ContainsKey("output_suffix"))
                input.Add("output_suffix", "");

            this.flexpde_script = (string)input["FlexPDE_file"];
            if (((string)input["output_suffix"]).Length == 4)
                this.flexpde_output = flexpde_script.Remove(flexpde_script.Length - 4) + "_old.pg6";
            else
            {
                string tmp_output =  (flexpde_script.Remove(flexpde_script.Length - 4) + (string)input["output_suffix"]);
                this.flexpde_output = tmp_output.Remove(tmp_output.Length - 4) + ".pg6";
            }

            // for flexpde, we need the -S parameter included in the input
            flexpde_options = "-S " + flexpde_script;

            // set the tolerances for FlexPDE if specified in the input dictionary
            if (input.ContainsKey("pot_tol"))
                pot_tol = (double)input["pot_tol"];
            if (input.ContainsKey("newton_tol"))
                newton_tol = (double)input["newton_tol"];
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
