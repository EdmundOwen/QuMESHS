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

namespace Solver_Bases.Mixing_Schedulers
{
    class File_Defined_Scheduler : Scheduler_Base
    {
        int[] switching_times;
        double[] alpha_vals;
        double[] zeta_vals;

        public File_Defined_Scheduler(Dictionary<string, object> dict)
        {
            // load data from file
            string[] data = File.ReadAllLines((string)dict["filename"]);
            // and find where the schedule data starts and finishes
            int start = -1; int end = -1;
            for (int i = 0; i < data.Length; i++)
            {
                if (data[i].Trim() == "{")
                    start = i + 1;
                if (data[i].Trim() == "}")
                    end = i - 1;
            }
            // check that the data has been found
            if (start == -1 || end == -1)
                throw new Exception("Error - cannot find schedule data in file: " + (string)dict["filename"]);

            // instantiate the arrays
            switching_times = new int[end - start + 1];
            alpha_vals = new double[end - start + 1];
            zeta_vals = new double[end - start + 1];

            // and load the data
            for (int i = 0; i < end - start + 1; i++)
            {
                string[] split_data = data[start + i].Split(' ');
                switching_times[i] = int.Parse(split_data[0].Split('=')[1]);
                alpha_vals[i] = double.Parse(split_data[1].Split('=')[1]);
                zeta_vals[i] = double.Parse(split_data[2].Split('=')[1]);
            }
        }

        public override double[] Get_Mixing_Parameter(int count)
        {
            for (int i = 0; i < switching_times.Length - 1; i++)
                if (count >= switching_times[i] && count < switching_times[i + 1])
                    return new double[] { alpha_vals[i], zeta_vals[i] };

            throw new IndexOutOfRangeException();
        }
    }
}
