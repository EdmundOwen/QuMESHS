using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases.Layers
{
    static public class Layer_Tool
    {
        static public string Get_Permitivity(Material material)
        {
            switch (material)
            {
                case Material.GaAs:
                    return "eps_r_GaAs * eps_0";
                case Material.AlGaAs:
                    return "eps_r_AlGaAs * eps_0";
                case Material.InGaAs:
                    return "eps_r_InGaAs * eps_0";
                case Material.InAlAs:
                    return "eps_r_InAlAs * eps_0";
                case Material.PMMA:
                    return "eps_pmma * eps_0";
                case Material.Metal:
                    return "0.0";
                case Material.Air:
                    return "eps_0";
                default:
                    throw new NotImplementedException("Error - Cannot find electrical permitivity for material: " + material.ToString());
            }
        }

        static public Material GetMaterial(string material)
        {
            switch (material)
            {
                case "air":
                    return Material.Air;

                case "metal":
                    return Material.Metal;

                case "gaas":
                    return Material.GaAs;

                case "algaas":
                    return Material.AlGaAs;

                case "ingaas":
                    return Material.InGaAs;

                case "inalas":
                    return Material.InAlAs;

                case "pmma":
                    return Material.PMMA;

                case "substrate":
                    return Material.Substrate;

                default:
                    throw new NotImplementedException("Error - Material properties not known");
            }
        }

        /// <summary>
        /// gets the layer numbers which have non-zero donor concentrations
        /// </summary>
        public static int[] Get_Donor_Layers(ILayer[] layers)
        {
            int[] result = (from item in layers
                               where item.Acceptor_Conc != 0.0 || item.Donor_Conc != 0.0
                               select item.Layer_No).ToArray();

            return result;
        }
    }
}
