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
