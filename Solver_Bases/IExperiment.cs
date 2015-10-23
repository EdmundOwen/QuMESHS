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
using Solver_Bases.Layers;

namespace Solver_Bases
{
    public interface IExperiment
    {
        bool Run();
        void Initialise(Dictionary<string, object> input_dict);

        SpinResolved_Data Carrier_Density { get; }
        SpinResolved_Data Dopent_Density { get; }
        Band_Data Chemical_Potential { get; }
        Band_Data GPhi { get; }
        Band_Data X { get; }
        ILayer[] Layers { get; }

        double Temperature { get; }
        double Current_Temperature { get; }

        int Nx_Dens { get; }
        double Dx_Dens { get; }
        double Xmin_Dens { get; }
        int Ny_Dens { get; }
        double Dy_Dens { get; }
        double Ymin_Dens { get; }
        int Nz_Dens { get; }
        double Dz_Dens { get; }
        double Zmin_Dens { get; }

        int Nx_Pot { get; }
        double Dx_Pot { get; }
        double Xmin_Pot { get; }
        int Ny_Pot { get; }
        double Dy_Pot { get; }
        double Ymin_Pot { get; }
        int Nz_Pot { get; }
        double Dz_Pot { get; }
        double Zmin_Pot { get; }
    }
}
