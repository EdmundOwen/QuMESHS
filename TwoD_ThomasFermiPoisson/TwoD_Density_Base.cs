using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases;
using Solver_Bases.Geometry;
using Solver_Bases.Layers;
using CenterSpace.NMath.Core;

namespace TwoD_ThomasFermiPoisson
{
    public abstract class TwoD_Density_Base : Density_Base
    {
        IExperiment exp;
        protected double dx, dy;
        protected double xmin, ymin;
        protected int nx, ny;

        protected Plane plane;
        protected double pos_z;

        public TwoD_Density_Base(IExperiment exp) : this(exp, Plane.yz, 0.0)
        {
        }

        public TwoD_Density_Base(IExperiment exp, Plane dir, double plane_pos)
            : base(exp.Temperature)
        {
            this.exp = exp;
            this.pos_z = plane_pos;

            // Get nx and ny dimensions depending on the orientation
            if (plane == Plane.xy)
            {
                this.dx = exp.Dx_Dens; this.dy = exp.Dy_Dens;
                this.xmin = exp.Xmin_Dens; this.ymin = exp.Ymin_Dens;
                this.nx = exp.Nx_Dens; this.ny = exp.Ny_Dens;
            }
            else if (plane == Plane.yz)
            {
                this.dx = exp.Dy_Dens; this.dy = exp.Dz_Dens;
                this.xmin = exp.Ymin_Dens; this.ymin = exp.Zmin_Dens;
                this.nx = exp.Ny_Dens; this.ny = exp.Nz_Dens;
            }
            else if (plane == Plane.zx)
            {
                this.dx = exp.Dz_Dens; this.dy = exp.Dx_Dens;
                this.xmin = exp.Zmin_Dens; this.ymin = exp.Xmin_Dens;
                this.nx = exp.Nz_Dens; this.ny = exp.Nx_Dens;
            }
            else throw new NotImplementedException();
        }

        public override SpinResolved_Data Get_ChargeDensity_Deriv(ILayer[] layers, SpinResolved_Data carrier_density_deriv, SpinResolved_Data dopent_density_deriv, Band_Data chem_pot)
        {
            // artificially deepen the copies of spin up and spin down
            Band_Data tmp_spinup = new Band_Data(new DoubleMatrix(carrier_density_deriv.Spin_Up.mat.Rows, carrier_density_deriv.Spin_Up.mat.Cols));
            Band_Data tmp_spindown = new Band_Data(new DoubleMatrix(carrier_density_deriv.Spin_Down.mat.Rows, carrier_density_deriv.Spin_Down.mat.Cols));

            for (int i = 0; i < carrier_density_deriv.Spin_Up.mat.Rows; i++)
                for (int j = 0; j < carrier_density_deriv.Spin_Up.mat.Cols; j++)
                {
                    tmp_spinup.mat[i, j] = carrier_density_deriv.Spin_Up.mat[i, j];
                    tmp_spindown.mat[i, j] = carrier_density_deriv.Spin_Down.mat[i, j];
                }

            SpinResolved_Data new_density = new SpinResolved_Data(tmp_spinup, tmp_spindown);

            // finally, get the charge density and send it to this new array
            Get_ChargeDensity_Deriv(layers, ref new_density, chem_pot);

            return new_density;
        }

        public override SpinResolved_Data Get_ChargeDensity(ILayer[] layers, SpinResolved_Data carrier_density, SpinResolved_Data dopent_density, Band_Data chem_pot)
        {
            // artificially deepen the copies of spin up and spin down
            Band_Data tmp_spinup = new Band_Data(new DoubleMatrix(carrier_density.Spin_Up.mat.Rows, carrier_density.Spin_Up.mat.Cols));
            Band_Data tmp_spindown = new Band_Data(new DoubleMatrix(carrier_density.Spin_Down.mat.Rows, carrier_density.Spin_Down.mat.Cols));

            for (int i = 0; i < carrier_density.Spin_Up.mat.Rows; i++)
                for (int j = 0; j < carrier_density.Spin_Up.mat.Cols; j++)
                {
                    tmp_spinup.mat[i, j] = carrier_density.Spin_Up.mat[i, j];
                    tmp_spindown.mat[i, j] = carrier_density.Spin_Down.mat[i, j];
                }

            SpinResolved_Data new_density = new SpinResolved_Data(tmp_spinup, tmp_spindown);

            // finally, get the charge density and send it to this new array
            Get_ChargeDensity(layers, ref new_density, chem_pot);

            return new_density;
        }

        public override double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input)
        {
            ZeroD_Density chem_pot_cal = new ZeroD_Density(Solver_Bases.Geometry.Geom_Tool.GetLayer(layers, z), temperature_input);

            return chem_pot_cal.Get_Equilibrium_Chemical_Potential();
        }

        public override void Close()
        {
            Console.WriteLine("Closing density solver");
        }

        protected void Get_Potential(ref Band_Data dft_band_offset, ILayer[] layers)
        {
            for (int i = 0; i < nx; i++)
                for (int j = 0; j < ny; j++)
                {
                    double pos_x = xmin + i * dx;
                    double pos_y = ymin + j * dy;
                    double band_gap = Geom_Tool.GetLayer(layers, plane, pos_x, pos_y, pos_z).Band_Gap;

                    dft_band_offset.mat[i, j] = 0.5 * band_gap - dft_band_offset.mat[i, j];
                }
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
        public abstract void Get_ChargeDensity_Deriv(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot);
    }
}
