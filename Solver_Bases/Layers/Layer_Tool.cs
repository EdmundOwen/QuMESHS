/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For additional information, please contact eto24@cam.ac.uk or visit
 * <http://www.qumeshs.org>
 * 
 **************************************************************************/

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
                case Material.Al2O3:
                    return "eps_al2o3 * eps_0";
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

                case "al2o3":
                    return Material.Al2O3;

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
