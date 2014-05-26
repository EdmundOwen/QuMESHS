using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;

namespace Solver_Bases.Geometry
{
    public static class Geom_Tool
    {
        public static ILayer GetLayer(ILayer[] layers, double z)
        {
            return GetLayer(layers, 0.0, 0.0, z);
        }

        public static ILayer GetLayer(ILayer[] layers, double y, double z)
        {
            return GetLayer(layers, 0.0, y, z);
        }

        public static ILayer GetLayer(ILayer[] layers, double x, double y, double z)
        {
            ILayer result = null;
            for (int i = 0; i < layers.Length; i++)
                if (layers[i].InLayer(x, y, z))
                    result = layers[i];
            
            // if the layer is not found, throw an exception
            if (result == null)
                throw new Exception("Error - layer not found at (x, y, z) = (" + x.ToString() + ", " + y.ToString() + ", " + z.ToString() + ")");

            // else if the result is the substrate and this is not the surface, throw an exception
            if (result.Material == Material.Substrate && z != result.Zmax)
                throw new Exception("Error - layer at (x, y, z) = (" + x.ToString() + ", " + y.ToString() + ", " + z.ToString() + ") is in the substrate");

            return result;
        }

        public static int Find_Layer_Below_Surface(ILayer[] layers)
        {
            for (int i = 0; i < layers.Length; i++)
                if (layers[i].Zmax == 0.0)
                    return layers[i].Layer_No;

            throw new Exception("Error - cannot find the layer immediately below the surface");
        }

        public static int Find_Layer_Above_Surface(ILayer[] layers)
        {
            for (int i = 0; i < layers.Length; i++)
                if (layers[i].Zmin == 0.0)
                    return layers[i].Layer_No;

            throw new Exception("Error - cannot find the layer immediately below the surface");
        }

        public static double Get_Xmin(ILayer[] layers)
        {
            double result = double.MaxValue;

            // ignore the substrate
            for (int i = 0; i < layers.Length; i++)
            {
                if (layers[i].Material == Material.Substrate)
                    continue;

                if (layers[i].Xmin < result)
                    result = layers[i].Xmin;
            }

            return result;
        }

        public static double Get_Ymin(ILayer[] layers)
        {
            double result = double.MaxValue;

            // only look through the top n-1 layers as you know the bottom will be the substrate
            for (int i = 0; i < layers.Length; i++)
            {
                if (layers[i].Material == Material.Substrate)
                    continue;

                if (layers[i].Ymin < result)
                    result = layers[i].Ymin;
            }

            return result;
        }

        public static double Get_Zmin(ILayer[] layers)
        {
            double result = double.MaxValue;

            // only look through the top n-1 layers as you know the bottom will be the substrate
            for (int i = 0; i < layers.Length; i++)
            {
                if (layers[i].Material == Material.Substrate)
                    continue;

                if (layers[i].Zmin < result)
                    result = layers[i].Zmin;
            }

            return result;
        }
    }
}
