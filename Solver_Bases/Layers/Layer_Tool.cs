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
                    return "eps_r * eps_0";
                case Material.AlGaAs:
                    return "eps_r * eps_0";
                case Material.PMMA:
                    return "eps_pmma * eps_0";
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
                case "gaas":
                    return Material.GaAs;

                case "algaas":
                    return Material.AlGaAs;

                case "pmma":
                    return Material.PMMA;

                case "substrate":
                    return Material.Substrate;

                default:
                    throw new NotImplementedException("Error - Material properties not known");
            }
        }
    }
}
