using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using Solver_Bases;
using Solver_Bases.Layers;

namespace Solver_GUI
{
    public partial class Solver_GUI : Form
    {
        int dimension = 2;
        Band_Data chem_pot;
        Band_Data x;
        SpinResolved_Data car_dens;
        SpinResolved_Data dop_dens;
        IExperiment exp;

        public Solver_GUI()
        {
            InitializeComponent();
            this.dimensionality.SelectedIndex = 1;
        }

        private void exitToolStripMenuItem_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        private void dimensionality_SelectedIndexChanged(object sender, EventArgs e)
        {
            this.Text = "Solver - " + this.dimensionality.SelectedItem.ToString();
            this.dimension = this.dimensionality.SelectedIndex + 1;
            
            // check dimensionality to alter inputs
            if (dimension == 1)
            {
                this.nx1Dval.Enabled = false;
                this.ny1Dval.Enabled = false;
            }
            else if (dimension == 2)
            {
                this.nx1Dval.Enabled = false;
                this.ny1Dval.Enabled = true;
            }
            else
            {
                this.nx1Dval.Enabled = true;
                this.ny1Dval.Enabled = true;
            }
        }

        private void dftCheck_CheckedChanged(object sender, EventArgs e)
        {
            if (!dftCheck.Checked)
            {
                this.dftnz1Dval.Enabled = false;
                this.dftzmin1Dval.Enabled = false;
            }
            else
            {
                this.dftnz1Dval.Enabled = true;
                this.dftzmin1Dval.Enabled = true;
            }
        }

        int i = 0;
        private void run_button_Click(object sender, EventArgs e)
        {
            conduction_band.Series["conduction_band_data"].Points.Clear();
            conduction_band.Series["x_data"].Points.Clear();
            density.Series["car_dens_data"].Points.Clear();
            density.Series["dop_dens_data"].Points.Clear();
            density.Series["gphi_data"].Points.Clear();
 //           density.Series["rho_prime_data"].Points.Clear();

            if (car_dens == null)
            {
                car_dens = new SpinResolved_Data(1900);
                dop_dens = new SpinResolved_Data(1900);
            }
       //     for (int i = 0; i < 200; i++)
       //    {
                exp = new OneD_ThomasFermiPoisson.Experiment();

                Dictionary<string, object> inputs = new Dictionary<string, object>();
                Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Solver_Config_tmp.txt");
                Inputs_to_Dictionary.Add_Input_Parameters_to_Dictionary(ref inputs, "Input_Parameters_1D_tmp.txt");
                inputs["surface_charge"] = 0.0;

                exp.Initialise(inputs);
                exp.Run();

                ILayer[] layers = exp.Layers;
                double dz = exp.Dz_Pot;
                int nz = exp.Nz_Pot;
                double zmin = exp.Zmin_Pot;
                x = exp.X;
                chem_pot = (Input_Band_Structure.Get_BandStructure_Grid(layers, dz, nz, zmin) - exp.Chemical_Potential + x);
                car_dens = exp.Carrier_Density;
                dop_dens = exp.Dopent_Density;
                Band_Data gphi = exp.GPhi;
            SpinResolved_Data rhop = exp.Rho_Prime;

                for (int j = 0; j < chem_pot.Length; j++)
                {
                    double pos = zmin + dz * j;
                    conduction_band.Series["conduction_band_data"].Points.AddXY(pos, chem_pot[j]);
                    conduction_band.Series["x_data"].Points.AddXY(pos, x[j]);
                    density.Series["car_dens_data"].Points.AddXY(pos, car_dens.Spin_Summed_Data[j]);
                    density.Series["dop_dens_data"].Points.AddXY(pos, dop_dens.Spin_Summed_Data[j]);
                    density.Series["gphi_data"].Points.AddXY(pos, gphi[j]);
      //              density.Series["rho_prime_data"].Points.AddXY(pos, rhop.Spin_Summed_Data[j]);
                }


                conduction_band.ChartAreas["ChartArea1"].AxisY.Minimum = -60.0;
                conduction_band.ChartAreas["ChartArea1"].AxisY.Maximum = 60.0;
                conduction_band.ChartAreas["ChartArea1"].AxisY.Interval = 20.0;

                conduction_band.Refresh();
                density.Refresh();
                this.count_label.Text = "Count = " + (i + 1).ToString();
                this.output_label.Text = exp.Output_Info.Replace("\t", "   ");
                i++;
                this.count_label.Refresh();
                exp.Initialise(inputs);
      //      }
            
        }

        private void refresh_button_Click(object sender, EventArgs e)
        {
            conduction_band.Series["conduction_band_data"].Points.Clear();
            conduction_band.Series["x_data"].Points.Clear();
            density.Series["car_dens_data"].Points.Clear();
            density.Series["dop_dens_data"].Points.Clear();
            density.Series["gphi_data"].Points.Clear();
       //     density.Series["rho_prime_data"].Points.Clear();
            i = 0;
            System.IO.File.Delete("restart.flag");
        }
    }
}
