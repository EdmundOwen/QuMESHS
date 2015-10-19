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
