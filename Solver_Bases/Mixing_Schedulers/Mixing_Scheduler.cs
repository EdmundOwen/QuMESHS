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
