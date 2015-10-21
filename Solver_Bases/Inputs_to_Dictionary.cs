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
