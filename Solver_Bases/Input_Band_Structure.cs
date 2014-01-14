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
    static public class Input_Band_Structure
    {
        /*
        double[] layer_boundaries;
        int end_of_band_structure;

        string[] raw_input;
        Material[] band_structure;

        string band_structure_filename;
        
        public Input_Band_Structure()
        {
            // get an array of the desired materials structure
            raw_input = Get_BandStructure_Data(filename, out layer_boundaries, out end_of_band_structure);
            band_structure = Get_Materials_Array(end_of_band_structure, raw_input);
        }
        */

        public static ILayer[] Get_Layers(string filename)
        {
            Dictionary<int, Dictionary<string, object>> data = Get_BandStructure_Data(filename);

            // find which layer is the surface (the first layer is the surface by default)
            int surface = -1; int no_layers = data.Count;
            for (int i = 0; i < data.Count; i++)
                if ((string)data[i]["raw_data"] == "surface=true")
                {
                    surface = i;
                    data[i].Add("is_surface", true);
                    no_layers--;
                }

            Unpack_RawData(data);
            Calculate_Layer_Positions(data, surface);

            ILayer[] layers = new ILayer[no_layers];
            int index = 0;
            for (int i = 0; i < data.Count; i++)
                if (i == surface)
                    continue;
                else
                {
                    layers[index] = Create_Layer(data[i]);
                    index++;
                }

            return layers;
        }

        static ILayer Create_Layer(Dictionary<string, object> data)
        {
            // get various necessary properties from this list
            Material mat = Layer_Tool.GetMaterial((string)data["mat"]);
            Geometry_Type geom_type;
            if (data.ContainsKey("geom")) { geom_type = GetGeometryType((string)data["geom"]); } else geom_type = Geometry_Type.slab;

            // set geometry
            IGeom geom;
            switch (geom_type)
            {
                case Geometry_Type.slab:
                    geom = new Slab((double)data["zmin"], (double)data["zmax"]);
                    break;

                case Geometry_Type.sheet:
                    geom = new Sheet((double)data["zmin"]);
                    break;

                default:
                    throw new NotImplementedException("Error - Unknown geometry");
            }

            // create layer
            ILayer result;
            switch (mat)
            {
                case Material.GaAs:
                    result = new GaAs_Layer(geom, (int)data["layer_no"]);
                    break;

                case Material.AlGaAs:
                    result = new AlGaAs_Layer(geom, (int)data["layer_no"]);
                    break;

                case Material.PMMA:
                    result = new PMMA_Layer(geom, (int)data["layer_no"]);
                    break;

                case Material.Air:
                    result = new Air_Layer(geom, (int)data["layer_no"]);
                    break;

                // substrate always comes from -infty
                case Material.Substrate:
                    result = new Substrate_Layer(new Slab(-1.0 * double.MaxValue, (double)data["zmax"]), (int)data["layer_no"]);
                    break;

                default:
                    throw new NotImplementedException("Error - Unknown material");
            }

            // and set dopent levels if necessary
            if (data.ContainsKey("nd") || data.ContainsKey("na"))
            {
                if (!data.ContainsKey("nd"))
                    data.Add("nd", 0.0);
                if (!data.ContainsKey("na"))
                    data.Add("na", 0.0);

                // use a factor of 10^-21 to convert from cm^-3 to nm^-3
                data["na"] = 1.0e-21 * (double)data["na"];
                data["nd"] = 1.0e-21 * (double)data["nd"];

                result.Set_Dopents((double)data["na"], (double)data["nd"]);
            }

            return result;
        }

        static void Unpack_RawData(Dictionary<int, Dictionary<string, object>> data)
        {
            for (int i = 0; i < data.Count; i++)
            {
                string[] raw_layer_data = ((string)data[i]["raw_data"]).Trim().Split(' ');
                for (int j = 0; j < raw_layer_data.Length; j++)
                {
                    string tmp_key = raw_layer_data[j].Split('=')[0].ToLower();
                    string tmp_value = raw_layer_data[j].Split('=')[1].ToLower();

                    // try and pass tmp_value by any means possible
                    double d_val; bool b_val;
                    if (double.TryParse(tmp_value, out d_val))
                        data[i].Add(tmp_key, d_val);
                    else if (bool.TryParse(tmp_value, out b_val))
                        data[i].Add(tmp_key, b_val);
                    else
                        data[i].Add(tmp_key, tmp_value);
                }

                data[i].Add("layer_no", i);
            }
        }

        static void Calculate_Layer_Positions(Dictionary<int, Dictionary<string, object>> data, int surface)
        {
            double tmp_position = 0.0;
            // calculate layer positions below the surface
            for (int i = surface - 1; i > 0; i--)
                if (!data[i].ContainsKey("zmin"))
                {
                    data[i].Add("zmax", tmp_position);
                    tmp_position -= (double)data[i]["t"];
                    data[i].Add("zmin", tmp_position);
                }
                else
                    tmp_position = (double)data[i]["zmin"];

            // and set the substrate
            data[0].Add("zmax", tmp_position);
            data[0].Add("zmin", -1.0 * double.MaxValue);

            tmp_position = 0.0;
            // and now above the surface
            for (int i = surface + 1; i < data.Count; i++)
                if (!data[i].ContainsKey("zmax"))
                {
                    data[i].Add("zmin", tmp_position);
                    tmp_position += (double)data[i]["t"];
                    data[i].Add("zmax", tmp_position);
                }
                else
                    tmp_position = (double)data[i]["zmax"];
        }

        static Geometry_Type GetGeometryType(string geometry_type)
        {
            switch (geometry_type)
            {
                case "slab":
                    return Geometry_Type.slab;

                case "sheet":
                    return Geometry_Type.sheet;

                default:
                    throw new NotImplementedException("Error - Geometry not known");
            }
        }

        static Dictionary<int, Dictionary<string, object>> Get_BandStructure_Data(string filename)
        {
            if (!File.Exists(filename))
                throw new FileNotFoundException();

            // read data from input file (discarding comment lines and white space)
            string[] raw_input = (from line in File.ReadAllLines(filename)
                                  where !line.StartsWith("#") && line.Trim().Length != 0
                                  select line).ToArray();

            // calculate where the keyword "substrate" is... this indicates the end of the band structure part of the file
            int end_of_band_structure = 0;
            for (int i = 0; i < raw_input.Length; i++)
                if (raw_input[i].Split(' ')[0].Trim() == "mat=substrate")
                {
                    end_of_band_structure = i + 1;
                    break;
                }

            // create a dictionary of dictionaries each containing data on the given layer
            Dictionary<int, Dictionary<string, object>> result = new Dictionary<int, Dictionary<string, object>>();

            // containing the trimmed data string and find out which layer is directly below the surface
            string[] layer_data = new string[end_of_band_structure];
            // note, the layer order is reversed here 
            for (int i = end_of_band_structure - 1; i >= 0; i--)
            {
                result[i] = new Dictionary<string, object>();
                result[i].Add("raw_data", raw_input[end_of_band_structure - 1 - i]);
            }

            return result;
        }

        /*
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
                if (concentration[i] != 0.0)
                    energy[i] = Get_Dopent_Energy(band_structure[layer - 1], dopent);
            }
        }
        
        /// <summary>
        /// converts a string array of inputs to an array of Materials enums
        /// </summary>
        Material[] Get_Materials_Array(int end_of_band_structure, string[] raw_input)
        {
            Material[] band_structure = new Material[end_of_band_structure - 1];
            for (int i = 1; i < end_of_band_structure; i++)
                band_structure[i - 1] = Get_Material(raw_input[i]);

            return band_structure;
        }

        /// <summary>
        /// gets a materials enum based a Snider-type string input
        /// </summary>
        Material Get_Material(string input_file_entry)
        {
            string[] material = input_file_entry.Split();

            switch (material[0])
            {
                case "GaAs":
                    return Material.GaAs;

                case "AlAs":
                    return Material.AlAs;

                case "AlGaAs":
                    for (int i = 0; i < material.Length; i++)
                        if (material[i].StartsWith("x"))
                            if (material[i].Split('=').Last() == ".3" || material[i].Split('=').Last() == "0.3")
                                return Material.Al03GaAs;
                    
                    goto default;

                case "InAs":
                    return Material.InAs;
                    
                case "InGaAs":
                    for (int i = 0; i < material.Length; i++)
                        if (material[i].StartsWith("x"))
                            if (material[i].Split('=').Last() == ".75" || material[i].Split('=').Last() == "0.75")
                                return Material.In075GaAs;
                    
                    goto default;

                case "PMMA":
                    return Material.PMMA;

                default:
                    throw new NotImplementedException("Error - Material properties not known for " + input_file_entry[0]);
            }
        }

        /// <summary>
        /// returns the band gap for the given material in meV
        /// </summary>
        public double Get_BandGap(Material materials)
        {
            switch (materials)
            {
                case Material.GaAs:
                    return 1420.0;
                case Material.Al03GaAs:
                    return 1800.0;
                case Material.PMMA:
                    return 4400.0;

                default:
                    throw new NotImplementedException("Error - no material details for " + materials.ToString());
            }
        }

        /// <summary>
        /// returns the dopent energy for the given material in meV above or below E_f
        /// </summary>
        public double Get_Dopent_Energy(Material materials, Dopent dopent)
        {
            if (dopent == Dopent.acceptor)
            {
                switch (materials)
                {
                    case Material.GaAs:
                        return 680.0;
                    case Material.Al03GaAs:
                        return 859.0;

                    default:
                        throw new NotImplementedException("Error - no material details for " + materials.ToString());
                }
            }
            else if (dopent == Dopent.donor)
            {
                switch (materials)
                {
                    case Material.GaAs:
                        return 704.0;
                    case Material.Al03GaAs:
                        return 869.0;

                    default:
                        throw new NotImplementedException("Error - no material details for " + materials.ToString());
                }
            }
            else throw new NotImplementedException();
        }
        */

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
        /// returns a DoubleMatrix with the given band structure planarised in the growth direction
        /// </summary>
        public static Band_Data Expand_BandStructure(DoubleVector structure, int nx, int ny)
        {
            DoubleMatrix[] result = new DoubleMatrix[structure.Length];

            for (int i = 0; i < structure.Length; i++)
            {
                result[i] = new DoubleMatrix(nx, ny);
                for (int j = 0; j < nx; j++)
                    for (int k = 0; k < ny; k++)
                        result[i][j, k] = structure[i];
            }

            return new Band_Data(result);
        }

        public static Band_Data Get_BandStructure_Grid(ILayer[] layers, double dz, int nz, double zmin)
        {
            DoubleVector result = new DoubleVector(nz);
            for (int i = 0; i < nz; i++)
                result[i] = 0.5 * Geom_Tool.GetLayer(layers, zmin + i * dz).Band_Gap;

            return new Band_Data(result);
        }

        public static Band_Data Get_BandStructure_Grid(ILayer[] layers, double dy, double dz, int ny, int nz, double ymin, double zmin)
        {
            DoubleMatrix result = new DoubleMatrix(ny, nz);
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                    result[i, j] = 0.5 * Geom_Tool.GetLayer(layers, ymin + i * dy, zmin + j * dz).Band_Gap;

            return new Band_Data(result);
        }

        public static Band_Data Get_BandStructure_Grid(ILayer[] layers, double dx, double dy, double dz, int nx, int ny, int nz, double xmin, double ymin, double zmin)
        {
            throw new NotImplementedException();
        }

        /*
        /// <summary>
        /// array of materials from shallowest to deepest
        /// </summary>
        public Material[] Materials_Array
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
        */
    }
}
