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
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public interface ILayer
    {
        bool Dopents_Frozen_Out(double temperature);
        void Set_Dopents(double acceptor_concentration, double donor_concentration);

        #region Access_Geometry
        bool InLayer(double z);
        bool InLayer(double y, double z);
        bool InLayer(double x, double y, double z);
        Geometry_Type Geometry { get; }
        ILayer Get_Component(int component_no);
        int No_Components { get; }

        double Xmin { get; }
        double Xmax { get; }
        double Ymin { get; }
        double Ymax { get; }
        double Zmin { get; }
        double Zmax { get; }
        #endregion

        #region Layer_Properties
        int Layer_No { get; }
        Material Material { get; }
        double Permitivity { get; }
        double Band_Gap { get; }
        double Donor_Energy { get; }
        double Acceptor_Energy { get; }
        double Donor_Conc { get; }
        double Acceptor_Conc { get; }
        double Dopent_FreezeOut_T { get; }
        #endregion

        double Electron_Mass { get; }
        double Hole_Mass { get; }
    }
}
