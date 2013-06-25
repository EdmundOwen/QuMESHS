using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Solver_Bases
{
    public abstract class Density_Solver : Physics_Base
    {
        protected double fermi_Energy, temperature;
        protected double dx, dy, dz;
        protected int nx, ny, nz;

        public Density_Solver(double fermi_Energy, double temperature, double dx, double dy, double dz, int nx, int ny, int nz)
        {
            this.fermi_Energy = fermi_Energy; this.temperature = temperature;
            this.dx = dx; this.dy = dy; this.dz = dz;
            this.nx = nx; this.ny = ny; this.nz = nz;
        }

        /// <summary>
        /// Gets the dopent occupation probability at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_NonSpinResolved_Dopent_Occupation(double energy)
        {
            return Get_NonSpinResolved_Dopent_Occupation(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Calculates the dopent occupation probability for arbitrary energy, E_f and T
        /// </summary>
        protected double Get_NonSpinResolved_Dopent_Occupation(double energy, double E_f, double T)
        {
            if (T == 0)
                if (energy > E_f)
                    return 0.0;
                else
                    return 1.0;
            else
                return 2.0 / (Math.Exp((energy - fermi_Energy) / (kB * T)) + 2.0);
        }

        /// <summary>
        /// Gets the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function(double energy)
        {
            return Get_Fermi_Function(energy, fermi_Energy, temperature);
        }
    }
}
