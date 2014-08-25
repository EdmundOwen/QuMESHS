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
        Experiment exp;

        public TwoD_Density_Base(Experiment exp)
            : base(exp.Temperature, 0.0, exp.Dy_Dens, exp.Dz_Dens, 1, exp.Ny_Dens, exp.Nz_Dens, 0.0, exp.Ymin_Dens, exp.Zmin_Dens)
        {
            this.exp = exp;
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
            for (int i = 0; i < ny; i++)
                for (int j = 0; j < nz; j++)
                {
                    double pos_y = ymin + i * dy;
                    double pos_z = zmin + j * dz;
                    double band_gap = Geom_Tool.GetLayer(layers, pos_y, pos_z).Band_Gap;
                    dft_band_offset.mat[i, j] = 0.5 * band_gap - dft_band_offset.mat[i, j];
                }
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
