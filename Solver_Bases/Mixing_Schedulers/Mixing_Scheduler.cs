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
    enum Scheduler_Type
    {
        None,
        Constant,
        File_Defined
    }

    public static class Mixing_Scheduler_Factory
    {
        public static IScheduler Get_Scheduler(string schedule_file)
        {
            // load the schedule file
            Dictionary<string, object> schedule_dict = new Dictionary<string,object>();
            if (File.Exists(schedule_file))
                Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref schedule_dict, schedule_file);
            else
                schedule_dict.Add("type", "none");

            // get the type of scheduler required and return the scheduler
            switch (((string)schedule_dict["type"]).ToLower())
            {
                case "none":
                    // these are the default values for the mixing scheduler
                    schedule_dict.Add("alpha", 0.01);
                    schedule_dict.Add("zeta", 3.0);
                    schedule_dict.Remove("type");
                    schedule_dict.Add("type", Scheduler_Type.None);
                    return new Constant_Scheduler(schedule_dict);

                case "constant":
                    schedule_dict.Remove("type");
                    schedule_dict.Add("type", Scheduler_Type.Constant);
                    return new Constant_Scheduler(schedule_dict);

                case "file_defined":
                    schedule_dict.Add("filename", schedule_file);
                    schedule_dict.Remove("type");
                    schedule_dict.Add("type", Scheduler_Type.File_Defined);
                    return new File_Defined_Scheduler(schedule_dict);

                default:
                    throw new NotImplementedException("Error - the scheduler type requested is not implemented");
            }
        }        
    }
}
