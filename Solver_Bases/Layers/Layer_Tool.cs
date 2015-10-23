/***************************************************************************
 * 
 * QuMESHS (Quantum Mesoscopic Electronic Semiconductor Heterostructure
 * Solver) for calculating electron and hole densities and electrostatic
 * potentials using self-consistent Poisson-Schroedinger solutions in 
 * layered semiconductors
 * 
 * Copyright(C) 2015 E. T. Owen and C. H. W. Barnes
 * 
 * The MIT License (MIT)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
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
