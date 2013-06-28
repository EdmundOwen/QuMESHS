using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public abstract class Potential_Base : Physics_Base
    {
        protected int nx, ny, nz;
        protected double dx, dy, dz;

        protected bool flexPDE;
        protected string flexpde_inputfile;

        public Potential_Base(double dx, double dy, double dz, int nx, int ny, int nz, bool using_flexPDE, string flexPDE_input)
        {
            this.nx = nx; this.ny = ny; this.nz = nz;
            this.dx = dx; this.dy = dy; this.dz = dz;

            // check whether using flexPDE
            flexPDE = using_flexPDE;
            if (using_flexPDE)
                this.flexpde_inputfile = flexPDE_input;
        }

        /// <summary>
        /// gets the potential using flexPDE
        /// </summary>
        protected Potential_Data Get_Potential_From_FlexPDE(Potential_Data density)
        {
            // save density to file in a FlexPDE "TABLE" format
            Save_Density(density, "density_1d.dat");

            // remove pot.dat if it still exists (to make sure that a new data file is made by flexPDE)
            try { File.Delete("pot.dat"); }
            catch (Exception) { }

            // run the flexPDE program as a process (quietly)
            Process.Start("C:\\FlexPDE6\\FlexPDE6.exe", "-Q " + flexpde_inputfile);
            while (!File.Exists("pot.dat"))
                Thread.Sleep(10);
            Thread.Sleep(10);

            string[] lines = File.ReadAllLines("pot.dat");

            // work out where the data starts (this is flexPDE specific)
            int first_line = 0;
            for (int i = 0; i < lines.Length; i++)
                if (lines[i].StartsWith("}"))
                {
                    first_line = i + 1;
                    break;
                }

            return Parse_Potential(lines, first_line);
        }

        /// <summary>
        /// Checks whether the potential has converged by comparing old and new potentials
        /// and determining whether every term is the same within a given tolerance
        /// </summary>
        public bool Check_Convergence(Potential_Data potential, Potential_Data new_potential, double tol)
        {
            Potential_Data pot_diff = new_potential - potential;

            int[] converged_test = new int[pot_diff.Length];
            for (int i = 0; i < nz; i++)
            {
                converged_test[i] = 0;
                if (Math.Abs(pot_diff[i]) < tol)
                    converged_test[i] = 1;
            }

            if (converged_test.Sum() == pot_diff.Length)
                return true;
            else
                return false;
        }

        abstract void Save_Density(Potential_Data density, string filename);
        abstract Potential_Data Parse_Potential(string[] data, int first_line);
    }
}
