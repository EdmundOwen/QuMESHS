using System;
using System.Collections.Generic;
using System.IO;
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

        public void Output(SpinResolved_DoubleVector data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            for (int i = 0; i < data.Nx; i++)
                sw.WriteLine(data.Spin_Summed_Vector[i].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_DoubleMatrix data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - Ordering compared to PotentialData objects is not guaranteed!");

            for (int i = 0; i < data.Nx; i++)
                for (int j = 0; j < data.Ny; j++)
                    sw.WriteLine(data.Spin_Summed_Matrix[i, j].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_DoubleMatrix[] data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - Ordering compared to PotentialData objects is not guaranteed!");

            for (int i = 0; i < data[0].Nx; i++)
                for (int j = 0; j < data[0].Ny; j++)
                    for (int k = 0; k < data.Length; k++)
                        sw.WriteLine(data[k].Spin_Summed_Matrix[i, j].ToString());

            sw.Close();
        }

        /// <summary>
        /// Gets the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function(double energy)
        {
            return Get_Fermi_Function(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Gets the fermi function for the dopents at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Dopent_Fermi_Function(double energy)
        {
            return Get_Dopent_Fermi_Function(energy, fermi_Energy, temperature);
        }
    }
}
