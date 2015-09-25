using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CenterSpace.NMath.Core;

namespace Solver_Bases
{
    public static class Physics_Base
    {
        // physical constants
        public const double hbar = 0.658211814;                 // (meV) (ps)
        public const double m_e = 5.68562958e-3;                // (meV) (ps)^2 (nm)^-2 of the free electron mass
        public const double q_e = 160.217646;                   // (zC) is positive as in these definitions it is the elementary charge
        public const double energy_V_to_meVpzC = 6.2415093;     // conversion factor from V to meV per zC
        public const double epsilon_0 = 1.41859713;             // (zC)^2 (nm)^-1 (meV)^-1 for vacuum
        public const double epsilon_r = 13.0;                   // relative permittivity for vacuum -> GaAs
        public const double epsilon_r_GaAs = 12.9;              // relative permittivity for vacuum -> GaAs
        public const double epsilon_r_InAs = 14.6;              // relative permittivity for vacuum -> InAs
        public const double epsilon_r_AlAs = 10.1;              // relative permittivity for vacuum -> AlAs
        public const double epsilon_pmma = 2.6;                 // relative permittivity for vacuum -> PMMA
        public const double epsilon = epsilon_0 * epsilon_r;    // (zC)^2 (nm)^-1 (meV)^-1 for GaAs
        public const double kB = 0.086173324;                   // (meV) (K)^-1

        // particle-specific constants (may be changed but only through Change_Mass(double new_mass) method
        // default values are for electrons in GaAs
        public static double a_B;       // Bohr radius for GaAs in nm
        public static double Ry;     // rydberg in meV in GaAs

        /// <summary>
        /// gets 3D spin-resolved density of states for a given potential and energy
        /// </summary>
        public static double Get_3D_DensityofStates(double energy, double band_edge, Carrier carrier_type, double mass)
        {
            // if the energy is below the potential (i.e. below the conduction band) return zero
            if (carrier_type == Carrier.Electron && energy < band_edge)
                return 0.0;
            else if (carrier_type == Carrier.Hole && energy > band_edge)
                return 0.0;

            // generate density of states
            if (carrier_type == Carrier.Electron)
                return Math.Pow(2.0 * mass * mass * mass * (energy - band_edge), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);
            else if (carrier_type == Carrier.Hole)
                return Math.Pow(2.0 * mass * mass * mass * (band_edge - energy), 0.5) / (Math.PI * Math.PI * hbar * hbar * hbar);
            else
                throw new NotImplementedException();
        }

        static double Get_Rs(double charge_density)
        {
            if (charge_density == 0.0)
                return double.PositiveInfinity;

            charge_density = Math.Abs(charge_density) / Physics_Base.q_e;    // DFT works on density, not charge 
            return Math.Pow(0.75 / (Math.PI * charge_density), 1.0 / 3.0) / Physics_Base.a_B;    // and convert to dimensionless constant
        }
/*
        static Dictionary<string, double> Unpolarized_Parameters = new Dictionary<string, double>
        {
            {"kappa", 0.9164},
            {"A", 0.0622},
            {"B", -0.096},
            {"gamma", -0.2846},
            {"beta1", 1.0529},
            {"beta2", 0.3334},
            {"C", 0.0040},
            {"D", -0.0232}
        };
        static Dictionary<string, double> Polarized_Parameters = new Dictionary<string, double>
        {
            {"kappa", 1.1546},
            {"A", 0.03110},
            {"B", -0.0538},
            {"gamma", -0.1686},
            {"beta1", 1.3981},
            {"beta2", 0.2611},
            {"C", 0.0014},
            {"D", -0.0096}
        };*/
        static Dictionary<string, double> Unpolarized_Parameters = new Dictionary<string, double>
        {
            {"A", 0.031091},
            {"B", 0.21370},
            {"C", 7.5957},
            {"D", 3.5876},
            {"E", 1.6382},
            {"F", 0.49294}
        };
        static Dictionary<string, double> Polarized_Parameters = new Dictionary<string, double>
        {
            {"A", 0.015545},
            {"B", 0.20548},
            {"C", 14.1189},
            {"D", 6.1977},
            {"E", 3.3662},
            {"F", 0.62517}
        };
        static Dictionary<string, double> Minus_Spin_Stiffness_Parameters = new Dictionary<string, double>
        {
            {"A", 0.016887},
            {"B", 0.11125},
            {"C", 10.357},
            {"D", 3.6231},
            {"E", 0.88026},
            {"F", 0.49671}
        };

        /// <summary>
        /// returns Perdew-Wang generalised function (see Guiliani & Vignale pp.49-50)
        /// </summary>
        static double Get_PW_G_Function(double r_s, Dictionary<string, double> p)
        {
            if (r_s == 0.0)
                return 0.0;
            else if (r_s < 0.0)
                throw new ArgumentOutOfRangeException("Error - r_s must be positive");
            else
            {
                double g = p["A"] * (1.0 + p["C"] * Math.Pow(r_s, 0.5) + p["D"] * r_s + p["E"] * Math.Pow(r_s, 1.5) + p["F"] * r_s * r_s);
                return -4.0 * p["A"] * (1.0 + p["B"] * r_s) * Math.Log(1.0 + 0.5 / g);
            }
        }
        static double Get_PW_G_Function_Deriv(double r_s, Dictionary<string, double> p)
        {
            if (r_s == 0.0)
                return 0.0;
            else if (r_s < 0.0)
                throw new ArgumentOutOfRangeException("Error - r_s must be positive");
            else
            {
                double g = p["A"] * (1.0 + p["C"] * Math.Pow(r_s, 0.5) + p["D"] * r_s + p["E"] * Math.Pow(r_s, 1.5) + p["F"] * r_s * r_s);
                double dg = p["A"] * (0.5 * p["C"] * Math.Pow(r_s, -0.5) + p["D"] + 1.5 * p["E"] * Math.Pow(r_s, 0.5) + 2.0 * p["F"] * r_s);

                return -4.0 * p["A"] * p["B"] * Math.Log(1.0 + 0.5 / g) 
                        - 4.0 * p["A"] * (1.0 + p["B"] * r_s) / (1.0 + 0.5 / g) * -0.5 * dg / (g * g);
            }
        }

        /// <summary>
        /// return the Perdew - Wang correlation potential
        /// </summary>
        public static double Get_LSDA_Ec_Potential(double dens_up, double dens_down)
        {
            double polarisation = (dens_up - dens_down) / (dens_up + dens_down);

            double ec, d_ec;                                                 // correlation potential
            double r_s = Get_Rs(dens_up + dens_down);
            
            if (double.IsInfinity(r_s))
                return 0.0;

            // correlation energy per particle from Guiliani & Vignale (Perdew-Wang form)
            ec = Get_PW_G_Function(r_s, Unpolarized_Parameters);
            if (polarisation != 0.0)
                ec += -1.0 * Get_PW_G_Function(r_s, Minus_Spin_Stiffness_Parameters) * Polarisation_Interpolation(polarisation) / 1.709921 * (1.0 - Math.Pow(polarisation, 4.0))
                        + (Get_PW_G_Function(r_s, Polarized_Parameters) - Get_PW_G_Function(r_s, Unpolarized_Parameters)) * Polarisation_Interpolation(polarisation) * Math.Pow(polarisation, 4.0);

            // and its derivative
            d_ec = Get_PW_G_Function_Deriv(r_s, Unpolarized_Parameters);
            if (polarisation != 0.0)
                d_ec += -1.0 * Get_PW_G_Function_Deriv(r_s, Minus_Spin_Stiffness_Parameters) * Polarisation_Interpolation(polarisation) / 1.709921 * (1.0 - Math.Pow(polarisation, 4.0))
                                + (Get_PW_G_Function_Deriv(r_s, Polarized_Parameters) - Get_PW_G_Function_Deriv(r_s, Unpolarized_Parameters)) * Polarisation_Interpolation(polarisation) * Math.Pow(polarisation, 4.0);

            // potential from LDA form
            double v_c = ec - r_s * d_ec / 3.0;

            // convert from Ry to meV
            v_c *= Physics_Base.Ry;

            return v_c;
        }

        static double Polarisation_Interpolation(double p)
        {
            return (Math.Pow(1.0 - p, 4.0 / 3.0) + Math.Pow(1.0 + p, 4.0 / 3.0) - 2.0) / (Math.Pow(2.0, 4.0 / 3.0) - 2.0);
        }

        public static double Get_LSDA_Ex_Potential(double dens_up, double dens_down)
        {
            double polarisation = (dens_up - dens_down) / (dens_up + dens_down);
            return 0.5 * Get_Ex_Potential(dens_up + dens_down, Unpolarized_Parameters) * (Math.Pow(1.0 + polarisation, 4.0 / 3.0) + Math.Pow(1.0 - polarisation, 4.0 / 3.0));
        }

        /// <summary>
        /// return the exchange potential (see Guiliani & Vignale p. 34)
        /// </summary>
        static double Get_Ex_Potential(double charge_density, Dictionary<string, double> xc_params)
        {
            double v_x;                                                 // exchange potential
            double r_s = Get_Rs(charge_density);

            if (double.IsInfinity(r_s))
                return 0.0;

            // exchange energy per particle for uniform electron gas from Guiliani & Vignale p. 34
            // e_x = - 0.9164 / r_s
            v_x = -1.0 * (4.0 * 0.9164 / 3.0) / r_s;

            // convert from Ry to meV
            v_x *= Physics_Base.Ry;

            return v_x;
        }



        /// <summary>
        /// return the derivative, with respect to the density, of the Perdew - Zunger exchange correlation potential
        /// </summary>
        public static double Get_XC_Potential_Deriv(double charge_density)
        {
            double d_vc, d_vx;                                                 // correlation and exchange potentials
            double r_s = Get_Rs(charge_density);

            double density = Math.Abs(charge_density) / Physics_Base.q_e;
            double d_rs = -1.0 * r_s / (3.0 * density); // derivative of rs with respect to the density

            if (double.IsInfinity(r_s))
                return 0.0;

            // exchange energy per particle for uniform electron gas from Parr and Yang (Density Functional Theory of Atoms and Molecules)
            // e_x = - 0.9164 / r_s
            d_vx = (4.0 * 0.9164 / 3.0) / (r_s * r_s);

            // correlation energy per particle from Perdew and Zunger (1981) App C, v_xc = (1 - r_s/3 d/dr_s) e_xc (in Ry, note that PZ-1981 is in Ha)
            // e_c = -0.2846 / (1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s)                           r_s > 1
            // e_c = 0.0622 * Math.Log(r_s) - 0.0960 - 0.0232 * r_s + 0.0040 * r_s * Math.Log(r_s)      r_s < 1

            if (r_s > 1)
                d_vc = -0.2846 * (1 + 7.0 * 1.0529 / 6.0 * (0.5 / Math.Sqrt(r_s)) + 4.0 * 0.3334 / 3.0) * Math.Pow(1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s, -2.0)
                                + -0.2846 * (1 + 7.0 * 1.0529 / 6.0 * Math.Sqrt(r_s) + 4.0 * 0.3334 / 3.0 * r_s) * -2.0 * (1.0529 * 0.5 / Math.Sqrt(r_s) + 0.3334) * Math.Pow(1.0 + 1.0529 * Math.Sqrt(r_s) + 0.3334 * r_s, -3.0);
            else
                d_vc = 0.0622 / r_s + (2.0 * 0.0040 / 3.0) * (1.0 + Math.Log(r_s)) + (2.0 * -0.0232 - 0.0040) / 3.0;

            // convert from Ry to meV
            d_vx *= Physics_Base.Ry;
            d_vc *= Physics_Base.Ry;

            return d_rs * (d_vx + d_vc);
        }

        public static double Get_LSDA_XC_Potential(double dens_up, double dens_down)
        {
            if (dens_up + dens_down == 0.0)
                return 0.0;
           
            return Get_LSDA_Ex_Potential(dens_up, dens_down) + Get_LSDA_Ec_Potential(dens_up, dens_down);
        }

        public static double Get_XC_Potential(double charge_density)
        {
            if (charge_density == 0.0)
                return 0.0;

            // by default, give the unpolarized exchange-correlation potential
            return Get_LSDA_Ex_Potential(0.5 * charge_density, 0.5 * charge_density) + Get_LSDA_Ec_Potential(0.5 * charge_density, 0.5 * charge_density);
        }

        public static Band_Data Get_XC_Potential(SpinResolved_Data charge_density)
        {
            Band_Data result;
            Band_Data charge_dens_spin_summed = charge_density.Spin_Summed_Data;
            int dim = charge_dens_spin_summed.Dimension;

            if (dim == 1)
            {
                int nx = charge_dens_spin_summed.vec.Length;
                result = new Band_Data(nx, 0.0);
                for (int i = 0; i < nx; i++)
                    result.vec[i] = Get_XC_Potential(charge_dens_spin_summed.vec[i]);

                return result;
            }
            else if (dim == 2)
            {
                int nx = charge_dens_spin_summed.mat.Rows;
                int ny = charge_dens_spin_summed.mat.Cols;
                result = new Band_Data(nx, ny, 0.0);
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result.mat[i, j] = Get_XC_Potential(charge_dens_spin_summed.mat[i, j]);

                return result;
            }
            else if (dim == 3)
            {
                int nx = charge_dens_spin_summed.vol[0].Rows;
                int ny = charge_dens_spin_summed.vol[0].Cols;
                int nz = charge_dens_spin_summed.vol.Length;
                result = new Band_Data(nx, ny, nz, 0.0);
                for (int k = 0; k < nz; k++)
                    for (int i = 0; i < nx; i++)
                        for (int j = 0; j < ny; j++)
                            result.vol[k][i, j] = Get_XC_Potential(charge_dens_spin_summed.vol[k][i, j]);

                return result;
            }
            else
                throw new NotImplementedException();
        }

        public static Band_Data Get_XC_Potential_Deriv(SpinResolved_Data charge_density)
        {
            Band_Data result;
            Band_Data charge_dens_spin_summed = charge_density.Spin_Summed_Data;
            int dim = charge_dens_spin_summed.Dimension;

            if (dim == 1)
            {
                int nx = charge_dens_spin_summed.vec.Length;
                result = new Band_Data(new DoubleVector(nx));
                for (int i = 0; i < nx; i++)
                    result.vec[i] = Get_XC_Potential_Deriv(charge_dens_spin_summed.vec[i]);

                return result;
            }
            else if (dim == 2)
            {
                int nx = charge_dens_spin_summed.mat.Rows;
                int ny = charge_dens_spin_summed.mat.Cols;
                result = new Band_Data(new DoubleMatrix(nx, ny));
                for (int i = 0; i < nx; i++)
                    for (int j = 0; j < ny; j++)
                        result.mat[i, j] = Get_XC_Potential_Deriv(charge_dens_spin_summed.mat[i, j]);

                return result;
            }
            else if (dim == 3)
            {
                throw new NotImplementedException();
            }
            else
                throw new NotImplementedException();
        }
    }
}
