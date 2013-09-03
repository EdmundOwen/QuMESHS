using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public class Input_Band_Structure
    {
        double[] layer_boundaries;
        int end_of_band_structure;

        string[] raw_input;
        Materials[] band_structure;

        public Input_Band_Structure(string filename)
        {
            // get an array of the desired materials structure
            raw_input = Get_BandStructure_Data(filename, out layer_boundaries, out end_of_band_structure);
            band_structure = Get_Materials_Array(end_of_band_structure, raw_input);
        }

        /// <summary>
        /// outputs a DoubleVector containing the given band structure
        /// </summary>
        /// <param name="nz">number of points in the growth direction for the output</param>
        /// <param name="dz">point separation for the output</param>
        /// <returns></returns>
        public DoubleVector GetBandStructure(int nz, double dz)
        {
            DoubleVector result = new DoubleVector(nz);

            // fill the result vector with band structures
            int layer = 0;
            for (int i = 0; i < nz; i++)
            {
                double z = i * dz;
                if (z > layer_boundaries[layer + 1]  && layer != band_structure.Length - 1)
                    layer++;

                result[i] = Get_BandGap(band_structure[layer]);
            }

            return result;
        }

        public void GetDopentData(int nz, double dz, Dopent dopent, out DoubleVector concentration, out DoubleVector energy)
        {
            concentration = new DoubleVector(nz);
            energy = new DoubleVector(nz);

            string dopent_string;
            if (dopent == Dopent.acceptor)
                dopent_string = "na";
            else if (dopent == Dopent.donor)
                dopent_string = "nd";
            else
                throw new Exception("It is completely impossible to get here");

            // get the dopent concentration
            int layer = 1;
            for (int i = 0; i < nz; i++)
            {
                double z = i * dz;
                if (z > layer_boundaries[layer] && layer != end_of_band_structure - 1)
                    layer++;

                double conc = (from data in raw_input[layer].Split()
                                       where data.ToLower().StartsWith(dopent_string)
                                       select double.Parse(data.Split('=').Last())).FirstOrDefault();

                // add concentration to result with a factor of 10^-21 to convert from cm^-3 to nm^-3
                concentration[i] = conc * 1e-21;
                // and get the dopent energy
                energy[i] = Get_Dopent_Energy(band_structure[layer - 1], dopent);
            }
        }

        string[] Get_BandStructure_Data(string filename, out double[] layer_boundaries , out int end_of_band_structure)
        {
            if (!File.Exists(filename))
                throw new FileNotFoundException();

            // read data from input file (discarding comment lines and white space)
            string[] raw_input = (from line in File.ReadAllLines(filename)
                                  where !line.StartsWith("#") && line.Trim().Length != 0
                                  select line).ToArray();

            // calculate where the keyword "substrate" is... this indicates the end of the band structure part of the file
            end_of_band_structure = 0;
            for (int i = 0; i < raw_input.Length; i++)
                if (raw_input[i].Trim() == "substrate")
                {
                    end_of_band_structure = i;
                    break;
                }

            // get an array of their given thicknesses
            layer_boundaries = new double[end_of_band_structure];
            layer_boundaries[0] = 0.0;
            for (int i = 1; i < end_of_band_structure; i++)
                layer_boundaries[i] = double.Parse(raw_input[i].Split()[1].Split('=').Last()) + layer_boundaries[i - 1];

            return raw_input;
        }

        /// <summary>
        /// converts a string array of inputs to an array of Materials enums
        /// </summary>
        Materials[] Get_Materials_Array(int end_of_band_structure, string[] raw_input)
        {
            Materials[] band_structure = new Materials[end_of_band_structure - 1];
            for (int i = 1; i < end_of_band_structure; i++)
                band_structure[i - 1] = Get_Material(raw_input[i]);

            return band_structure;
        }

        /// <summary>
        /// gets a materials enum based a Snider-type string input
        /// </summary>
        Materials Get_Material(string input_file_entry)
        {
            string[] material = input_file_entry.Split();

            switch (material[0])
            {
                case "GaAs":
                    return Materials.GaAs;

                case "AlAs":
                    return Materials.AlAs;

                case "AlGaAs":
                    for (int i = 0; i < material.Length; i++)
                        if (material[i].StartsWith("x"))
                            if (material[i].Split('=').Last() == ".3" || material[i].Split('=').Last() == "0.3")
                                return Materials.Al03GaAs;
                    
                    goto default;

                case "InAs":
                    return Materials.InAs;
                    
                case "InGaAs":
                    for (int i = 0; i < material.Length; i++)
                        if (material[i].StartsWith("x"))
                            if (material[i].Split('=').Last() == ".75" || material[i].Split('=').Last() == "0.75")
                                return Materials.In075GaAs;
                    
                    goto default;

                default:
                    throw new NotImplementedException("Error - Material properties not known for " + input_file_entry[0]);
            }
        }

        /// <summary>
        /// returns the band gap for the given material in meV
        /// </summary>
        public double Get_BandGap(Materials materials)
        {
            switch (materials)
            {
                case Materials.GaAs:
                    return 1420.0;
                case Materials.Al03GaAs:
                    return 1800.0;

                default:
                    throw new NotImplementedException("Error - no material details for " + materials.ToString());
            }
        }

        /// <summary>
        /// returns the dopent energy for the given material in meV above or below E_f
        /// </summary>
        public double Get_Dopent_Energy(Materials materials, Dopent dopent)
        {
            if (dopent == Dopent.acceptor)
            {
                switch (materials)
                {
                    case Materials.GaAs:
                        return 680.0;
                    case Materials.Al03GaAs:
                        return 859.0;

                    default:
                        throw new NotImplementedException("Error - no material details for " + materials.ToString());
                }
            }
            else if (dopent == Dopent.donor)
            {
                switch (materials)
                {
                    case Materials.GaAs:
                        return 704.0;
                    case Materials.Al03GaAs:
                        return 869.0;

                    default:
                        throw new NotImplementedException("Error - no material details for " + materials.ToString());
                }
            }
            else throw new NotImplementedException();
        }

        /// <summary>
        /// returns a DoubleMatrix with the given band structure planarised in the transverse direction
        /// </summary>
        public static Band_Data Expand_BandStructure(DoubleVector structure, int ny)
        {
            DoubleMatrix result = new DoubleMatrix(ny, structure.Length);
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < structure.Length; j++)
                    result[i, j] = structure[j];

            return new Band_Data(result);
        }

        /// <summary>
        /// returns a DoubleMatrix with the given band structure planarised in the transverse direction
        /// </summary>
        public static Band_Data Expand_BandStructure(DoubleVector structure, int nx, int ny)
        {
            DoubleMatrix[] result = new DoubleMatrix[nx];

            for (int i = 0; i < nx; i++)
            {
                result[i] = new DoubleMatrix(ny, structure.Length);
                for (int j = 0; j < ny; j++)
                    for (int k = 0; k < structure.Length; k++)
                        result[i][j, k] = structure[k];
            }

            return new Band_Data(result);
        }

        /// <summary>
        /// array of materials from shallowest to deepest
        /// </summary>
        public Materials[] Materials_Array
        {
            get { return band_structure; }
        }

        /// <summary>
        /// array of layer boundary depths
        /// </summary>
        public double[] Layer_Depths
        {
            get { return layer_boundaries; }
        }
    }
}
