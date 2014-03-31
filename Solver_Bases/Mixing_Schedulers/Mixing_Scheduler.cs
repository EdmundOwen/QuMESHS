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
