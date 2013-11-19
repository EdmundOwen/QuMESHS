using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;
using Solver_Bases.Layers;
using Solver_Bases.Geometry;

namespace Solver_Bases
{
    public abstract class Experiment_Base
    {
        protected SpinResolved_Data charge_density;
        protected Band_Data band_offset;
        protected ILayer[] layers;

        protected double dx, dy, dz;
        protected double xmin, ymin, zmin = -1.0;
        protected int nx, ny, nz;

        protected double alpha, alpha_prime, tol;

        protected bool using_flexPDE = false;
        protected string flexPDE_input;
        protected string flexPDE_location;

        protected double initial_temperature = 300.0;
        protected double temperature;

        public void Initialise(Dictionary<string, object> input_dict)
        {
            // solver inputs
            if (input_dict.ContainsKey("tolerance")) this.tol = (double)input_dict["tolerance"]; else throw new KeyNotFoundException("No solution tolerance in input dictionary!");
            if (input_dict.ContainsKey("alpha")) { this.alpha = (double)input_dict["alpha"]; alpha_prime = alpha; } else throw new KeyNotFoundException("No potential mixing parameter, alpha, in input dictionary!");
        
            // will not use FlexPDE unless told to
            if (input_dict.ContainsKey("use_FlexPDE")) this.using_flexPDE = bool.Parse((string)input_dict["use_FlexPDE"]); else using_flexPDE = false;
            // default input file for FlexPDE is called "default.pde"
            if (input_dict.ContainsKey("FlexPDE_file")) this.flexPDE_input = (string)input_dict["FlexPDE_file"]; else this.flexPDE_input = "default.pde";
            if (using_flexPDE) 
            { 
                // make sure that FlexPDE does exist at the specified location
                try { this.flexPDE_location = (string)input_dict["FlexPDE_location"]; }
                catch (Exception) { throw new Exception("Error - no location for FlexPDE executable supplied"); }
                if (!File.Exists(flexPDE_location))
                    throw new Exception("Error - FlexPDE executable file does not exist at location" + flexPDE_location + "!");
            }
            
            // physical inputs
            if (input_dict.ContainsKey("init_T")) this.initial_temperature = (double)input_dict["init_T"];
            if (input_dict.ContainsKey("T")) this.temperature = (double)input_dict["T"]; else throw new KeyNotFoundException("No temperature in input dictionary!");

            //// get the band structure
            if (input_dict.ContainsKey("BandStructure_File"))
                layers = Input_Band_Structure.Get_Layers((string)input_dict["BandStructure_File"]);
            else throw new KeyNotFoundException("No band structure file found in input dictionary!");

            // and find the domain minimum coordinate values
            xmin = Geom_Tool.Get_Xmin(layers);
            ymin = Geom_Tool.Get_Ymin(layers);
            zmin = Geom_Tool.Get_Zmin(layers);
        }

        public abstract void Run();

        /// <summary>
        /// returns a list of dopent freeze-out temperatures between the initial and final temperatures of the experiment
        /// and with the final temperature at the end
        /// </summary>
        protected double[] Freeze_Out_Temperatures()
        {
            double[] raw_temps = new double[layers.Length];
            for (int i = 0; i < layers.Length; i++)
                raw_temps[i] = layers[i].Dopent_FreezeOut_T;

            // now, calculate which temperatures are in the range (final_T < T < init_T) and sort in descending order
            double[] temp_list = (from value in raw_temps where (value > temperature && value < initial_temperature) select value).ToArray().Concat(new[] {temperature}).ToArray();
            return temp_list.Distinct().OrderByDescending(c => c).ToArray();
        }
    }
}
