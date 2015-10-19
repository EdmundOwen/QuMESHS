/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public static class Inputs_to_Dictionary
    {
        /// <summary>
        /// Adds the inputs from filename into the dictionary... this includes file names for the wavefunctions and potentials
        /// <\summary>
        public static void Add_Input_Parameters_to_Dictionary(ref Dictionary<string, object> inputs, string filename)
        {
            // check that the input file exists
            if (!File.Exists(filename))
                throw new FileNotFoundException("Error - cannot find the file " + filename + " at this location...");

            // read all data from input file
            string[] raw_input = (from line in File.ReadAllLines(filename)
                                  where line.StartsWith("%")
                                  select line).ToArray();

            // input entries into the dictionary from the raw data
            for (int i = 0; i < raw_input.Length; i++)
            {
                // trim %
                raw_input[i] = raw_input[i].TrimStart('%');

                // find the index of the "="
                int index = raw_input[i].IndexOf('=');

                // separate the line into values before and after the '='
                string key = raw_input[i].Substring(0, index - 1).Trim();
                string val_string = raw_input[i].Substring(index + 1).Trim();
                double val;
                
                // if the val_string can be cast as a double, do so... else, just pass through the string
                if (double.TryParse(val_string, out val))
                    inputs.Add(key, val);
                // or try to convert it into a bool
                else if (val_string.ToLower() == "true" || val_string.ToLower() == "false")
                    inputs.Add(key, bool.Parse(val_string));
                else
                    inputs.Add(key, val_string);
            }

            Dictionary<string, double> output = new Dictionary<string, double>();
            Parse_Voltages(inputs, output);
            foreach (KeyValuePair<string, double> val in output)
                inputs.Add(val.Key, val.Value);
        }

        public static void Parse_Voltages(Dictionary<string, object> input, Dictionary<string, double> output)
        {
            // create a list of voltages (if they're available)
            if (input.ContainsKey("voltages"))
            {
                int count = 0;
                string[] raw_voltages = ((string)input["voltages"]).TrimStart('{').TrimEnd('}').Split(',');
                for (int i = 0; i < raw_voltages.Length; i++)
                {
                    output["V" + count.ToString()] = double.Parse(raw_voltages[i]);
                    count++;
                }
            }
        }

        public static void Output_Dictionary_To_File(Dictionary<string, object> input, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            foreach (KeyValuePair<string, object> option in input)
            {
                if (option.Value is bool || option.Value is int || option.Value is double || option.Value is string)
                    sw.WriteLine("% " + option.Key + " = " + option.Value.ToString());
            }
        }
    }
}
