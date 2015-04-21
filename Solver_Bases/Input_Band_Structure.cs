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
    public static class Input_Band_Structure
    {
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

                case Geometry_Type.strip:
                    geom = new Strip((double)data["zmin"], (double)data["zmax"], (double)data["dx"], (double)data["dy"], (double)data["width"], (double)data["theta"]);
                    break;

                case Geometry_Type.half_slab:
                    geom = new Half_Slab((double)data["zmin"], (double)data["zmax"], (double)data["dx"], (double)data["dy"], (double)data["theta"]);
                    break;

                case Geometry_Type.triangle_slab:
                    geom = new Triangle_Slab((double)data["zmin"], (double)data["zmax"], (double)data["x0"], (double)data["y0"], (double)data["theta1"], (double)data["theta2"]);
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
                    result = new AlGaAs_Layer(geom, (int)data["layer_no"], (double)data["x"]);
                    break;

                case Material.InGaAs:
                    result = new InGaAs_Layer(geom, (int)data["layer_no"], (double)data["x"]);
                    break;

                case Material.InAlAs:
                    result = new InAlAs_Layer(geom, (int)data["layer_no"], (double)data["x"]);
                    break;

                case Material.PMMA:
                    result = new PMMA_Layer(geom, (int)data["layer_no"]);
                    break;

                case Material.Metal:
                    result = new Metal_Layer(geom, (int)data["layer_no"]);
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

            // if this is actually a composite layer, then rewrite result... (this is simpler than putting a "Composite" material for the switch above)
            if (data.ContainsKey("composite"))
                if ((bool)data["composite"])
                    result = Create_Composite_Layer(result, data, geom);

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

        private static ILayer Create_Composite_Layer(ILayer default_layer, Dictionary<string, object> data, IGeom default_geom)
        {
            int no_components = (int)(double)data["no_components"];
            ILayer[] result = new ILayer[no_components];

            // set the default layer
            result[0] = default_layer;

            // unpack the rest of the composite data
            Dictionary<int, Dictionary<string, object>> composite_data = new Dictionary<int, Dictionary<string, object>>();
            Unpack_CompositeData(composite_data, data, default_geom);

            // and create the new layers from this
            for (int i = 1; i < no_components; i++)
                result[i] = Create_Layer(composite_data[i]);

            return new Composite_Layer(result, (int)data["layer_no"]);
        }

        static void Unpack_CompositeData(Dictionary<int, Dictionary<string, object>> composite_data, Dictionary<string, object> data, IGeom default_geom)
        {
            // cycle over the composite data (minus the default layer)
            for (int i = 1; i < (int)(double)data["no_components"]; i++)
            {
                string[] raw_component_data = ((string)data["component" + i]).TrimStart('{').TrimEnd('}').Split(',');

                Dictionary<string, object> component_data = new Dictionary<string, object>();

                // get the component data
                for (int j = 0; j < raw_component_data.Length; j++)
                {
                    string tmp_key = raw_component_data[j].Split('=')[0].ToLower();
                    string tmp_value = raw_component_data[j].Split('=')[1].ToLower();

                    // try and pass tmp_value by any means possible
                    double d_val; bool b_val;
                    if (double.TryParse(tmp_value, out d_val))
                        component_data.Add(tmp_key, d_val);
                    else if (bool.TryParse(tmp_value, out b_val))
                        component_data.Add(tmp_key, b_val);
                    else
                        component_data.Add(tmp_key, tmp_value);
                }

                // and add some geometry data
                component_data.Add("zmin", default_geom.Zmin);
                component_data.Add("zmax", default_geom.Zmax);
                component_data.Add("layer_no", i);

                composite_data[i] = component_data;
            }
        }

        static void Unpack_RawData(Dictionary<int, Dictionary<string, object>> data)
        {
            int layer_no = 1;

            for (int i = 0; i < data.Count; i++)
            {
                string[] raw_layer_data = ((string)data[i]["raw_data"]).Trim().Split(' ');
                int component_no = 1;

                for (int j = 0; j < raw_layer_data.Length; j++)
                {
                    //check if this is a composite layer, if so... just input the data in its raw format (ie. don't break it up into key/value pairs)
                    if (raw_layer_data[j].StartsWith("{"))
                    {
                        data[i].Add("component" + component_no, raw_layer_data[j]);
                        component_no++;
                        continue;
                    }

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

                data[i].Add("layer_no", layer_no);

                layer_no++;
                if (data[i].ContainsKey("is_surface"))
                    if ((bool)data[i]["is_surface"] == true)
                    {
                        data[i].Remove("layer_no");
                        layer_no--;
                    }
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

                case "strip":
                    return Geometry_Type.strip;

                case "half":
                    return Geometry_Type.half_slab;

                case "tria":
                    return Geometry_Type.triangle_slab;

                default:
                    throw new NotImplementedException("Error - Geometry not known");
            }
        }

        static Dictionary<int, Dictionary<string, object>> Get_BandStructure_Data(string filename)
        {
            if (!File.Exists(filename))
                throw new FileNotFoundException("Error - Could not find file \"" + filename + "\"");

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

        public static SpinResolved_Data Expand_BandStructure(SpinResolved_Data data, int ny)
        {
            Band_Data spin_up = Expand_BandStructure(data.Spin_Up.vec, ny);
            Band_Data spin_down = Expand_BandStructure(data.Spin_Down.vec, ny);

            return new SpinResolved_Data(spin_up, spin_down);
        }

        public static SpinResolved_Data Expand_BandStructure(SpinResolved_Data data, int nx, int ny)
        {
            Band_Data spin_up = Expand_BandStructure(data.Spin_Up.vec, nx, ny);
            Band_Data spin_down = Expand_BandStructure(data.Spin_Down.vec, nx, ny);

            return new SpinResolved_Data(spin_up, spin_down);
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
            DoubleMatrix[] result = new DoubleMatrix[nz];

            for (int k = 0; k < nz; k++)
            {
                result[k] = new DoubleMatrix(nx, ny);
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result[k][i, j] = 0.5 * Geom_Tool.GetLayer(layers, xmin + i * dx, ymin + j * dy, zmin + k * dz).Band_Gap;
            }

            return new Band_Data(result);
        }
    }
}
