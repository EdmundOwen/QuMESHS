using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

namespace Iterative_Greens_Function
{
    static class Inputs_to_Dictionary
    {
        /// <summary>
        /// Adds the inputs from filename into the dictionary... this includes file names for the wavefunctions and potentials
        /// <\summary>
        public static void Add_Input_Parameters_to_Dictionary(ref Dictionary<string, object> inputs, string filename)
        {
            inputs = new Dictionary<string, object>();

            // check that the input file exists
            if (!File.Exists(filename))
                throw new FileNotFoundException();

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
                else
                    inputs.Add(key, val_string);
            }
        }
    }
}
