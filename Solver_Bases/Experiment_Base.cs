using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public abstract class Experiment_Base
    {
        protected double alpha, alpha_prime, tol;

        protected bool using_flexPDE = false;
        protected string flexPDE_input;
        protected string flexPDE_location;

        protected Band_Data band_offset;

        protected double[] layer_depths;
        protected DoubleVector band_structure;
        protected DoubleVector acceptor_energy, donor_energy;
        protected DoubleVector acceptor_conc, donor_conc;

        protected double temperature = 300.0;

        // temperature at which the dopents are assumed to have frozen out
        protected double freeze_out_T = 70.0;
        protected bool freeze_dopents;

        public void Initialise(Dictionary<string, object> input_dict, double dz, int nz)
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
            if (input_dict.ContainsKey("T")) this.temperature = (double)input_dict["T"]; else throw new KeyNotFoundException("No temperature in input dictionary!");
            // and check whether the dopents are frozen out
            if (input_dict.ContainsKey("dopents_frozen")) this.freeze_dopents = (bool)input_dict["dopents_frozen"];
            else { if (temperature < freeze_out_T) freeze_dopents = true; else freeze_dopents = false; }

            //// get the band structure
            if (input_dict.ContainsKey("BandStructure_File"))
            {
                Input_Band_Structure band_structure_generator = new Input_Band_Structure((string)input_dict["BandStructure_File"]);
            
                band_structure = band_structure_generator.GetBandStructure(nz, dz);
                layer_depths = band_structure_generator.Layer_Depths;
            
                band_structure_generator.GetDopentData(nz, dz, Dopent.acceptor, out acceptor_conc, out acceptor_energy);
                band_structure_generator.GetDopentData(nz, dz, Dopent.donor, out donor_conc, out donor_energy);
            
                // finally, check that the total layer depth from the band structure generator is the same as nz
                if (layer_depths[layer_depths.Length - 1] != nz * dz)
                    throw new Exception("Error - The band structure specified is not the same as the given sample depth!");
            }
            else throw new KeyNotFoundException("No band structure file found in input dictionary!");
        }

        public abstract void Run();
    }
}
