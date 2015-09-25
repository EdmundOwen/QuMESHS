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
using System.IO;
using CenterSpace.NMath.Core;

namespace Solver_GUI
{
    public partial class Solver_GUI : Form
    {
        int dimension = 2;
        Band_Data chem_pot;
        Band_Data val_pot;
        SpinResolved_Data car_dens;
        SpinResolved_Data dop_dens;
        IExperiment exp1d;
        Dictionary<string, object> inputs = new Dictionary<string,object>();

        public Solver_GUI()
        {
            InitializeComponent();
            this.dimensionality.SelectedIndex = 1;
            Update_BandStructure();
            this.bandstructureCombo.SelectedIndex = 0;
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
            if (!dft1DCheck.Checked)
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

        private void step_button_Click(object sender, EventArgs e)
        {
            dimension = 1;

            conduction_band.Series["conduction_band_data"].Points.Clear();
            conduction_band.Series["valence_band_data"].Points.Clear();
            conduction_band.Series["x_data"].Points.Clear();
            density.Series["car_dens_data"].Points.Clear();
            density.Series["dop_dens_data"].Points.Clear();
            density.Series["gphi_data"].Points.Clear();

            int nz;
            int i = int.Parse(count_no_label.Text);

            if (!int.TryParse(nz1Dval.Text, out nz))
            {
                MessageBox.Show("Format of Nz for dopent potential is invalid");
                return;
            }

            // initialise dopent experiment if the count is zero
            if (i == 0)
            {
                exp1d = new OneD_ThomasFermiPoisson.Experiment();

                car_dens = new SpinResolved_Data(nz);
                dop_dens = new SpinResolved_Data(nz);

                Dictionary<string, object> inputs_tmp = Create_Dictionary();
                foreach (KeyValuePair<string, object> option in inputs_tmp) inputs.Add(option.Key.Replace("_1d", ""), option.Value);
                

                inputs["surface_charge"] = 0.0;
                inputs["max_iterations"] = 0.0;

                exp1d.Initialise(inputs);
            }

            converged = exp1d.Run();

            ILayer[] layers = exp1d.Layers;
            double dz = exp1d.Dz_Pot;
            double zmin = exp1d.Zmin_Pot;

            car_dens = exp1d.Carrier_Density;
            dop_dens = exp1d.Dopent_Density;
            Band_Data x = Physics_Base.q_e * exp1d.X;
            Band_Data gphi = exp1d.GPhi;
            chem_pot = (Input_Band_Structure.Get_BandStructure_Grid(layers, dz, nz, zmin) - exp1d.Chemical_Potential);// + x);
            val_pot = (-1.0 * Input_Band_Structure.Get_BandStructure_Grid(layers, dz, nz, zmin) - exp1d.Chemical_Potential);// + x);

            for (int j = 0; j < chem_pot.Length; j++)
            {
                double pos = zmin + dz * j;
                conduction_band.Series["conduction_band_data"].Points.AddXY(pos, chem_pot[j]);
                conduction_band.Series["valence_band_data"].Points.AddXY(pos, val_pot[j]);
                conduction_band.Series["x_data"].Points.AddXY(pos, x[j]);
                density.Series["car_dens_data"].Points.AddXY(pos, car_dens.Spin_Summed_Data[j]);
                density.Series["dop_dens_data"].Points.AddXY(pos, dop_dens.Spin_Summed_Data[j]);
                density.Series["gphi_data"].Points.AddXY(pos, gphi[j]);
            }

            Set_Plot_Axes(dens_xmin_val.Text, dens_xmax_val.Text, density.ChartAreas["ChartArea1"].AxisX);
            Set_Plot_Axes(dens_ymin_val.Text, dens_ymax_val.Text, density.ChartAreas["ChartArea1"].AxisY);
            Set_Plot_Axes(pot_xmin_val.Text, pot_xmax_val.Text, conduction_band.ChartAreas["ChartArea1"].AxisX);
            Set_Plot_Axes(pot_ymin_val.Text, pot_ymax_val.Text, conduction_band.ChartAreas["ChartArea1"].AxisY);

            conduction_band.Refresh();
            density.Refresh();
            this.count_no_label.Text = (i + 1).ToString();
            this.temperature_val_label.Text = exp1d.Current_Temperature.ToString() + " K";

            this.carrier_dopent_density_Text.Text = (from val in car_dens.Spin_Summed_Data.vec
                                                     where val < 0.0
                                                     select -1.0e14 * val * dz / Physics_Base.q_e).ToArray().Sum().ToString("e3");
        }

        private void Set_Plot_Axes(string min, string max, System.Windows.Forms.DataVisualization.Charting.Axis axis)
        {
            double min_val, max_val;

            if (min != "")
                if (!double.TryParse(min, out min_val))
                {
                    MessageBox.Show("Error - " + min + " is not a valid format for axis plotting");
                }
                else
                    axis.Minimum = min_val;
            else
                axis.Minimum = double.NaN;

            if (max != "")
                if (!double.TryParse(max, out max_val))
                {
                    MessageBox.Show("Error - " + max + " is not a valid format for axis plotting");
                }
                else
                    axis.Maximum = max_val;
            else
                axis.Maximum = double.NaN;
        }

        private void refresh_button_Click(object sender, EventArgs e)
        {
            inputs = new Dictionary<string, object>();

            conduction_band.Series["conduction_band_data"].Points.Clear();
            conduction_band.Series["valence_band_data"].Points.Clear();
            conduction_band.Series["x_data"].Points.Clear();
            density.Series["car_dens_data"].Points.Clear();
            density.Series["dop_dens_data"].Points.Clear();
            density.Series["gphi_data"].Points.Clear();
            System.IO.File.Delete("restart.flag");

            count_no_label.Text = "0";
            temperature_val_label.Text = "N/A";
        }

        bool converged = false;
        private void run_button_Click(object sender, EventArgs e)
        {
            converged = false;
            while (!converged)
            {
                step_button_Click(sender, e);
                count_no_label.Refresh();
                temperature_val_label.Refresh();
            }
        }

        private Dictionary<string, object> Create_Dictionary()
        {
            Dictionary<string, object> result = new Dictionary<string, object>();

            result.Add("nz_1d", double.Parse(nz1Dval.Text));
            result.Add("dz_1d", double.Parse(dz1Dval.Text));

            result.Add("TF_only_1d", !dft1DCheck.Checked);
            result.Add("nz_dens_1d", double.Parse(dftnz1Dval.Text));
            result.Add("dz_dens_1d", (double)result["dz_1d"]);
            result.Add("zmin_dens_1d", double.Parse(dftzmin1Dval.Text));

            result.Add("T_1d", double.Parse(temperature1Dval.Text));
            result.Add("top_V_1d", double.Parse(topV1Dval.Text));
            result.Add("bottom_V_1d", double.Parse(bottomV1Dval.Text));

            result.Add("tolerance_1d", double.Parse(tolerance1Dval.Text));
            result.Add("max_iterations_1d", double.Parse(max_iterations1d_val.Text));
            result.Add("use_FlexPDE_1d", false);

            result.Add("BandStructure_File", bandstructurefilename.Text);

            return result;
        }

        private void addlayer_Button_Click(object sender, EventArgs e)
        {
            if (material_combo.SelectedItem == null)
                return;

            ListViewItem new_layer = new ListViewItem(new string[] {
                                                        (string)material_combo.SelectedItem,
                                                        newlayer_thicknessval.Text,
                                                        newlayer_xval.Text,
                                                        newlayer_ndval.Text,
                                                        newlayer_naval.Text}, -1);

            // insert layer into listview
            bandstructure_list.Items.Insert(bandstructureCombo.SelectedIndex + 1, new_layer); 
            Update_BandStructure();
        }

        private void Update_BandStructure()
        {
            double thickness = 0.0;
            // update bandstructure combo list
            bandstructureCombo.Items.Clear();
            for (int i = 0; i < bandstructure_list.Items.Count; i++)
            {
                bandstructureCombo.Items.Add(i.ToString() + " " + bandstructure_list.Items[i].Text);
                if (bandstructure_list.Items[i].SubItems[1].Text != "")
                    thickness += double.Parse(bandstructure_list.Items[i].SubItems[1].Text);
            }

            total_thickness_val.Text = thickness.ToString();
            dz1Dval.Text = (thickness / double.Parse(nz1Dval.Text)).ToString();
        }
    
        private void material_combo_SelectedIndexChanged(object sender, EventArgs e)
        {
            if ((string)material_combo.SelectedItem == "AlGaAs" || (string)material_combo.SelectedItem == "InGaAs" || (string)material_combo.SelectedItem == "InAlAs")
            {
                newlayer_xval.Enabled = true;
                newlayer_xval.Text = "0.33";
            }
            else
            {
                newlayer_xval.Enabled = false;
                newlayer_xval.Text = "";
            }
        }

        private void deleteLayer_Button_Click(object sender, EventArgs e)
        {
            if ((string)bandstructureCombo.SelectedItem == "surface" || (string)bandstructureCombo.SelectedItem == "substrate")
                return;

            bandstructure_list.Items.RemoveAt(bandstructureCombo.SelectedIndex);
            Update_BandStructure();
        }

        private void editLayer_Button_Click(object sender, EventArgs e)
        {
            if (material_combo.SelectedItem == null || bandstructureCombo.SelectedIndex == -1)
                return;

            ListViewItem new_layer = new ListViewItem(new string[] {
                                                        (string)material_combo.SelectedItem,
                                                        newlayer_thicknessval.Text,
                                                        newlayer_xval.Text,
                                                        newlayer_ndval.Text,
                                                        newlayer_naval.Text}, -1);

            // insert layer into listview
            bandstructure_list.Items.Insert(bandstructureCombo.SelectedIndex + 1, new_layer);
            bandstructure_list.Items.RemoveAt(bandstructureCombo.SelectedIndex);
            Update_BandStructure();
        }

        private void bandstructure_create_button_Click(object sender, EventArgs e)
        {
            StreamWriter sw = new StreamWriter(bandstructurefilename.Text);

            // create band structure file
            for (int i = 0; i < bandstructure_list.Items.Count; i++)
            {
                string layerstring;
                ListViewItem item = bandstructure_list.Items[i];

                if (bandstructure_list.Items[i].Text == "surface")
                {
                    sw.WriteLine("surface=true");
                    continue;
                }

                // add material property to layer
                layerstring = "mat=" + item.SubItems[0].Text;
                if (item.SubItems[0].Text != "substrate" && item.SubItems[1].Text == "")
                {
                    MessageBox.Show("Error - thickness of layer " + i.ToString() + " is undefined");
                    sw.Close();
                    File.Delete(bandstructurefilename.Text);
                    return;
                }
                else if (item.SubItems[0].Text != "substrate")
                // add thickness to layer
                layerstring += " t=" + item.SubItems[1].Text;

                // add alloy composition to layer
                if (item.SubItems[2].Text != "")
                    layerstring += " x=" + item.SubItems[2].Text;

                // add donors to layer
                if (item.SubItems[3].Text != "")
                    layerstring += " Nd=" + item.SubItems[3].Text;

                // add acceptors to layer
                if (item.SubItems[4].Text != "")
                    layerstring += " Na=" + item.SubItems[4].Text;

                sw.WriteLine(layerstring);
            }

            sw.Close();
        }

        private void mnuOpenBandStructure_Click(object sender, EventArgs e)
        {
            // load the band data from a given file
            openFileDialog1.InitialDirectory = Environment.CurrentDirectory;
            openFileDialog1.FileName = bandstructurefilename.Text;
            openFileDialog1.ShowDialog();
            // read data from input file (discarding comment lines and white space)
            string[] bandstructure = (from line in File.ReadAllLines(openFileDialog1.FileName)
                                  where !line.StartsWith("#") && line.Trim().Length != 0
                                  select line).ToArray();


            bandstructure_list.Items.Clear();
            for (int i = 0; i < bandstructure.Length; i++)
            {
                // get the layer data and put it into a dictionary
                string[] layer_data = bandstructure[i].Split();
                Dictionary<string, string> layer_dict = new Dictionary<string, string>();
                for (int j = 0; j < layer_data.Length; j++)
                    layer_dict.Add(layer_data[j].Split('=')[0], layer_data[j].Split('=')[1]);

                string mat = "";
                string t = "";
                string x = "";
                string nd = "";
                string na = "";
                ListViewItem layer;

                // check if it is the surface
                if (layer_dict.ContainsKey("surface"))
                {
                    layer = new ListViewItem(new System.Windows.Forms.ListViewItem.ListViewSubItem[] {
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "surface", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Window, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0))))}, -1);
                    layer.UseItemStyleForSubItems = false;
                }
                else if (layer_dict["mat"] == "substrate")
                {
                    if (layer_dict.ContainsKey("Na"))
                        na = layer_dict["Na"];
                    if (layer_dict.ContainsKey("Nd"))
                        nd = layer_dict["Nd"];

                    layer = new ListViewItem(new System.Windows.Forms.ListViewItem.ListViewSubItem[] {
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "substrate"),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, nd),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, na)}, -1);
                    layer.UseItemStyleForSubItems = false;
                }
                else
                {
                    mat = layer_dict["mat"];
                    t = layer_dict["t"];
                    if (layer_dict.ContainsKey("x"))
                        x = layer_dict["x"];
                    if (layer_dict.ContainsKey("Na"))
                        na = layer_dict["Na"];
                    if (layer_dict.ContainsKey("Nd"))
                        nd = layer_dict["Nd"];
                 
                    // or the substrate
                    layer = new ListViewItem(new string[] { mat, t, x, nd, na });
                }

                bandstructure_list.Items.Add(layer);
            }

            bandstructure_list.Refresh();
            Update_BandStructure();
        }
    }
}
