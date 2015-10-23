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

namespace Solver_Bases
{
    public enum Spin
    {
        Up = 0,
        Down = 1
    }

    public enum Carrier
    {
        electron,
        hole,
        both
    }

    public enum Direction
    {
        x = 0,
        y = 1,
        z = 2
    }

    public enum Plane
    {
        xy = 0,
        yz = 1,
        zx = 2
    }

    public enum Material
    {
        GaAs,
        Al03GaAs,
        AlGaAs,
        AlAs,
        In075GaAs,
        InAlAs,
        InGaAs,
        InAs,
        PMMA,
        Al2O3,
        Air,
        Metal,
        Substrate
    }

    public enum Dopent
    {
        donor,
        acceptor
    }

    public enum Geometry_Type
    {
        triangle_slab,
        half_slab,
        strip,
        half_strip,
        slab,
        sheet, 
        composite
    }

    public enum OneD_Density
    {
        thomasfermi,
        dft
    }

    public enum TwoD_Density
    {
        effectiveband,
        thomasfermi,
        sodft,
        dft
    }

    public enum ThreeD_Density
    {
        effectiveband,
        thomasfermi,
        twodthomasfermi_oneddft,
        iterativegreensfunction
    }
}
