using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    class Composite_Layer : ILayer
    {
        int layer_no;
        ILayer[] layers;
        ILayer default_layer;

        public Composite_Layer(ILayer[] layers_in, int layer_no)
        {
            this.layers = layers_in;
            default_layer = layers[0];

            // the default layer is at index 0 in layers_in and must be of Geometry.Slab type
            if (default_layer.Geometry != Geometry_Type.slab)
                throw new FormatException("Error - the default input geometry must be a slab!\nLayer at index 0 is of type " + layers_in[0].Geometry.ToString());

            // this is a layer, so check that all input layers have the same zmin and zmax
            for (int i = 0; i < layers.Length; i++)
            {
                if (default_layer.Zmin != layers[i].Zmin)
                    throw new Exception("Error - Composite layer does not have consistent bottom");
                if (default_layer.Zmax != layers[i].Zmax)
                    throw new Exception("Error - Composite layer does not have consistent top");
            }

            // and assign layer number and limits in the z direction
            this.layer_no = layer_no;
        }

        /// <summary>
        /// returns the default layer (ie. the one first defined in the composite layer input, or marked "default")
        /// </summary>
        public ILayer GetLayer(double z)
        {
            return default_layer;
        }

        public ILayer GetLayer(double y, double z)
        {
            return GetLayer(0.0, y, z);
        }

        /// <summary>
        /// returns the layer at (x, y, z), which will be the default if there are no other layers at this point
        /// should check for consistency to make sure two layers are not defined here
        /// </summary>
        public ILayer GetLayer(double x, double y, double z)
        {
            List<ILayer> layer_list = GetLayerList(layers, x, y, z);

            if (layer_list.Count > 2)
                throw new Exception("Error - more than two layers found at (x, y, z) = (" + x.ToString() + ", " + y.ToString() + ", " + z.ToString() + ")\nCheck consistency in band structure input");
            else if (layer_list.Count == 2)
            {
                if (layer_list[0] == default_layer)
                    return layer_list[1];
                else if (layer_list[1] == default_layer)
                    return layer_list[0];
                else throw new Exception("Error - more than two layers found at (x, y, z) = (" + x.ToString() + ", " + y.ToString() + ", " + z.ToString() + ")\nCheck consistency in band structure input");
            }
            else if (layer_list.Count == 1)
                return default_layer;
            else
                throw new Exception("Error - Layer not found at at (x, y, z) = (" + x.ToString() + ", " + y.ToString() + ", " + z.ToString());
        }

        public List<ILayer> GetLayerList(ILayer[] layers, double x, double y, double z)
        {
            List<ILayer> layers_found = new List<ILayer>();

            for (int i = 0; i < layers.Length; i++)
                if (layers[i].InLayer(x, y, z))
                    layers_found.Add(layers[i]);

            return layers_found;
        }

        /// <summary>
        /// returns the component... this should be in the order in the band structure input
        /// default layer is 0
        /// </summary>
        public ILayer Get_Component(int component_no)
        {
            return layers[component_no];
        }

        public int Layer_No
        {
            get { return layer_no; }
        }

        public int No_Components
        {
            get { return layers.Length; }
        }

        public bool Dopents_Frozen_Out(double temperature)
        {
            return default_layer.Dopents_Frozen_Out(temperature);
        }

        public void Set_Dopents(double acceptor_concentration, double donor_concentration)
        {
            this.default_layer.Set_Dopents(acceptor_concentration, donor_concentration);
        }

        public bool InLayer(double z)
        {
            return default_layer.InLayer(z);
        }

        public bool InLayer(double y, double z)
        {
            return default_layer.InLayer(y, z);
        }

        public bool InLayer(double x, double y, double z)
        {
            return default_layer.InLayer(x, y, z);
        }

        public Geometry_Type Geometry
        {
            get { return Geometry_Type.composite; }
        }

        public double Xmin
        {
            get { return default_layer.Xmin; }
        }

        public double Xmax
        {
            get { return default_layer.Xmax; }
        }

        public double Ymin
        {
            get { return default_layer.Ymin; }
        }

        public double Ymax
        {
            get { return default_layer.Ymax; }
        }

        public double Zmin
        {
            get { return default_layer.Zmin; }
        }

        public double Zmax
        {
            get { return default_layer.Zmax; }
        }

        public Material Material
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the material"); }
        }

        public double Permitivity
        {
            get { return default_layer.Permitivity; }
        }

        public double Band_Gap
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the band gap"); }
        }

        public double Donor_Energy
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the donor energy"); }
        }

        public double Acceptor_Energy
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the acceptor energy"); }
        }

        public double Donor_Conc
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the donor concentration"); }
        }

        public double Acceptor_Conc
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the acceptor concentration"); }
        }

        public double Dopent_FreezeOut_T
        {
            get { throw new FormatException("Error - This is a composite layer! Use the Get_Layer(x, y, z) methods to find the dopent freeze-out temperature"); }
        }
    }
}
