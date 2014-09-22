using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;

namespace ThreeD_SchrodingerPoissonSolver
{
    public abstract class ThreeD_Density_Base : Density_Base
    {
        IExperiment exp;
        protected double dx, dy, dz;
        protected double xmin, ymin, zmin;
        protected int nx, ny, nz;

        public ThreeD_Density_Base(IExperiment exp)
            : base(exp.Temperature)
        {
            this.exp = exp;

            this.dx = exp.Dx_Dens; this.dy = exp.Dy_Dens; this.dz = exp.Dz_Dens;
            this.xmin = exp.Xmin_Dens; this.ymin = exp.Ymin_Dens; this.zmin = exp.Zmin_Dens;
            this.nx = exp.Nx_Dens; this.ny = exp.Ny_Dens; this.nz = exp.Nz_Dens;
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            throw new NotImplementedException();
        }

        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            // artificially deepen the copies of spin up and spin down
            Band_Data tmp_spinup = new Band_Data(carrier_density.Spin_Up.vol[0].Rows, carrier_density.Spin_Up.vol[0].Cols, carrier_density.Spin_Up.vol.Length, 0.0);
            Band_Data tmp_spindown = new Band_Data(carrier_density.Spin_Down.vol[0].Rows, carrier_density.Spin_Down.vol[0].Cols, carrier_density.Spin_Down.vol.Length, 0.0);

            for (int k = 0; k < carrier_density.Spin_Up.vol.Length; k++)
            {
                // fill with data
                for (int i = 0; i < carrier_density.Spin_Up.vol[0].Rows; i++)
                    for (int j = 0; j < carrier_density.Spin_Up.vol[0].Cols; j++)
                    {
                        tmp_spinup.vol[k][i, j] = carrier_density.Spin_Up.vol[k][i, j];
                        tmp_spindown.vol[k][i, j] = carrier_density.Spin_Down.vol[k][i, j];
                    }
            }

            SpinResolved_Data new_density = new SpinResolved_Data(tmp_spinup, tmp_spindown);

            // finally, get the charge density and send it to this new array
            Get_ChargeDensity(layers, ref new_density, chem_pot);

            return new_density;
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            // artificially deepen the copies of spin up and spin down
            Band_Data tmp_spinup = new Band_Data(carrier_density_deriv.Spin_Up.vol[0].Rows, carrier_density_deriv.Spin_Up.vol[0].Cols, carrier_density_deriv.Spin_Up.vol.Length, 0.0);
            Band_Data tmp_spindown = new Band_Data(carrier_density_deriv.Spin_Down.vol[0].Rows, carrier_density_deriv.Spin_Down.vol[0].Cols, carrier_density_deriv.Spin_Down.vol.Length, 0.0);

            for (int k = 0; k < carrier_density_deriv.Spin_Up.vol.Length; k++)
            {
                // fill with data
                for (int i = 0; i < carrier_density_deriv.Spin_Up.vol[0].Rows; i++)
                    for (int j = 0; j < carrier_density_deriv.Spin_Up.vol[0].Cols; j++)
                    {
                        tmp_spinup.vol[k][i, j] = carrier_density_deriv.Spin_Up.vol[k][i, j];
                        tmp_spindown.vol[k][i, j] = carrier_density_deriv.Spin_Down.vol[k][i, j];
                    }
            }

            SpinResolved_Data new_density = new SpinResolved_Data(tmp_spinup, tmp_spindown);

            // finally, get the charge density and send it to this new array
            Get_ChargeDensity_Deriv(layers, ref new_density, chem_pot);

            return new_density;
        }

        public override double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, x, y, z), temperature_input);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }

        public override void Close()
        {
            Console.WriteLine("Closing density solver");
        }

        double xmin_pot;
        public double Xmin_Pot
        {
            set { xmin_pot = value; }
        }
        double dx_pot;
        public double Dx_Pot
        {
            set { dx_pot = value; }
        }
        double ymin_pot;
        public double Ymin_Pot
        {
            set { ymin_pot = value; }
        }
        double dy_pot;
        public double Dy_Pot
        {
            set { dy_pot = value; }
        }
        double zmin_pot;
        public double Zmin_Pot
        {
            set { zmin_pot = value; }
        }
        double dz_pot;
        public double Dz_Pot
        {
            set { dz_pot = value; }
        }
        public abstract void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
    }
}
