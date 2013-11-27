using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;

namespace Solver_Bases
{
    public abstract class Density_Base : IDensity_Solve
    {
        protected double temperature;

        protected double dx, dy, dz;
        protected double xmin, ymin, zmin;
        protected int nx, ny, nz;

        protected double fermi_Energy = 0.0;
        protected int dim;

        public Density_Base(double temperature, double dx, double dy, double dz, int nx, int ny, int nz, double xmin, double ymin, double zmin)
        {
            this.temperature = temperature;
            this.dx = dx; this.dy = dy; this.dz = dz;
            this.nx = nx; this.ny = ny; this.nz = nz;
            this.xmin = xmin; this.ymin = ymin; this.zmin = zmin;

            // check how many actual dimensions are present
            if (nx != 1)
                dim = 3;
            else if (ny != 1)
                dim = 2;
            else
                dim = 1;
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
                return 2.0 / (Math.Exp((energy - fermi_Energy) / (Physics_Base.kB * T)) + 2.0);
        }

        /*public void Output(SpinResolved_DoubleVector data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);

            // output the charge density
            for (int i = 0; i < data.Nx; i++)
                sw.WriteLine(data.Spin_Summed_Vector[i].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_DoubleMatrix data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - Ordering compared to Band_Data objects is not guaranteed!");
            
            // output the charge density
            for (int i = 0; i < data.Nx; i++)
                for (int j = 0; j < data.Ny; j++)
                    sw.WriteLine(data.Spin_Summed_Matrix[i, j].ToString());

            sw.Close();
        }

        public void Output(SpinResolved_DoubleMatrix[] data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - Ordering compared to Band_Data objects is not guaranteed!");

            // output the charge density
            for (int i = 0; i < data[0].Nx; i++)
                for (int j = 0; j < data[0].Ny; j++)
                    for (int k = 0; k < data.Length; k++)
                        sw.WriteLine(data[k].Spin_Summed_Matrix[i, j].ToString());

            sw.Close();
        }*/

        public void Output(SpinResolved_Data data, string filename)
        {
            StreamWriter sw = new StreamWriter(filename);
            sw.WriteLine("Warning - The data has been written out serially and there is no information as to which order the dimensions come in.");
            sw.WriteLine("Warning - Ordering compared to Band_Data objects is not guaranteed!");

            // output the charge density
            Band_Data tot_charge = data.Spin_Summed_Data;
            for (int i = 0; i < tot_charge.Length; i++)
                sw.WriteLine(tot_charge[i].ToString());

            sw.Close();
        }

        /// <summary>
        /// Gets the fermi function at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Fermi_Function(double energy)
        {
            return Physics_Base.Get_Fermi_Function(energy, fermi_Energy, temperature);
        }

        /// <summary>
        /// Gets the fermi function for the dopents at a given energy using the default temperature and fermi energy
        /// </summary>
        protected double Get_Dopent_Fermi_Function(double energy)
        {
            return Physics_Base.Get_Dopent_Fermi_Function(energy, fermi_Energy, temperature);
        }

        public virtual double Get_Chemical_Potential(double z, ILayer[] layers)
        {
            return Get_Chemical_Potential(z, layers, temperature);
        }

        public virtual double Get_Chemical_Potential(double y, double z, ILayer[] layers)
        {
            return Get_Chemical_Potential(y, z, layers, temperature);
        }

        public virtual double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers)
        {
            return Get_Chemical_Potential(x, y, z, layers, temperature);
        }

        /// <summary>
        /// set the default chemical potential to be calculated at (0,0,z).  Everything else requires this method to be overridden
        /// </summary>
        public virtual double Get_Chemical_Potential(double z, ILayer[] layers, double temperature_input)
        {
            return Get_Chemical_Potential(0.0, 0.0, z, layers, temperature_input);
        }

        /// <summary>
        /// set the default chemical potential to be calculated at (0,y,z).  Everything else requires this method to be overridden
        /// </summary>
        public virtual double Get_Chemical_Potential(double y, double z, ILayer[] layers, double temperature_input)
        {
            return Get_Chemical_Potential(0.0, y, z, layers, temperature_input);
        }

        public abstract void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
        public abstract double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input);
        public abstract void Close();
    }
}
