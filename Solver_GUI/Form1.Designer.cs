using System.Drawing;

namespace Solver_GUI
{
    partial class Solver_GUI
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            System.Windows.Forms.ListViewItem listViewItem1 = new System.Windows.Forms.ListViewItem(new System.Windows.Forms.ListViewItem.ListViewSubItem[] {
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "surface", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Window, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0))))}, -1);
            System.Windows.Forms.ListViewItem listViewItem2 = new System.Windows.Forms.ListViewItem(new System.Windows.Forms.ListViewItem.ListViewSubItem[] {
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "substrate"),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "", System.Drawing.SystemColors.WindowText, System.Drawing.SystemColors.Control, new System.Drawing.Font("Microsoft Sans Serif", 8.25F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)))),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "0.0"),
            new System.Windows.Forms.ListViewItem.ListViewSubItem(null, "0.0")}, -1);
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea1 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend1 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series1 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series2 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series3 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea2 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Legend legend2 = new System.Windows.Forms.DataVisualization.Charting.Legend();
            System.Windows.Forms.DataVisualization.Charting.Series series4 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series5 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series6 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.menuStrip1 = new System.Windows.Forms.MenuStrip();
            this.fileToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuNew = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuNewBandStructure = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuNewInputFile = new System.Windows.Forms.ToolStripMenuItem();
            this.solverConfigToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuOpen = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuOpenBandStructure = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuOpenInputFile = new System.Windows.Forms.ToolStripMenuItem();
            this.solverConfigToolStripMenuItem1 = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem1 = new System.Windows.Forms.ToolStripSeparator();
            this.mnuSave = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuSaveAs = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripMenuItem2 = new System.Windows.Forms.ToolStripSeparator();
            this.mnuExit = new System.Windows.Forms.ToolStripMenuItem();
            this.viewToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.settingToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.helpToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.gettingStartedToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.aboutToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.splitContainer = new System.Windows.Forms.SplitContainer();
            this.inputTab = new System.Windows.Forms.TabControl();
            this.bandstructure = new System.Windows.Forms.TabPage();
            this.total_thickness_label = new System.Windows.Forms.Label();
            this.total_thickness_val = new System.Windows.Forms.TextBox();
            this.bandstructurefilename = new System.Windows.Forms.TextBox();
            this.bandstructurefile_label = new System.Windows.Forms.Label();
            this.bandstructure_create_button = new System.Windows.Forms.Button();
            this.group_newlayer = new System.Windows.Forms.GroupBox();
            this.newlayer_naval = new System.Windows.Forms.TextBox();
            this.newlayer_na_label = new System.Windows.Forms.Label();
            this.newlayer_ndval = new System.Windows.Forms.TextBox();
            this.newlayer_nd_label = new System.Windows.Forms.Label();
            this.newlayer_xval = new System.Windows.Forms.TextBox();
            this.newlayer_x_label = new System.Windows.Forms.Label();
            this.newlayer_thicknessval = new System.Windows.Forms.TextBox();
            this.newlayer_thickness_label = new System.Windows.Forms.Label();
            this.material_label = new System.Windows.Forms.Label();
            this.material_combo = new System.Windows.Forms.ComboBox();
            this.layerselect_label = new System.Windows.Forms.Label();
            this.bandstructure_list = new System.Windows.Forms.ListView();
            this.material_header = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.thickness_header = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.alloy_header = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.donor_header = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.acceptor_header = ((System.Windows.Forms.ColumnHeader)(new System.Windows.Forms.ColumnHeader()));
            this.editLayer_Button = new System.Windows.Forms.Button();
            this.deleteLayer_Button = new System.Windows.Forms.Button();
            this.addlayer_Button = new System.Windows.Forms.Button();
            this.bandstructureCombo_label = new System.Windows.Forms.Label();
            this.bandstructureCombo = new System.Windows.Forms.ComboBox();
            this.potentialinputs = new System.Windows.Forms.TabPage();
            this.densityinputs = new System.Windows.Forms.TabPage();
            this.output_label = new System.Windows.Forms.Label();
            this.dopentinputs = new System.Windows.Forms.TabPage();
            this.groupPlotting = new System.Windows.Forms.GroupBox();
            this.pot_ymax_label = new System.Windows.Forms.Label();
            this.pot_ymax_val = new System.Windows.Forms.TextBox();
            this.pot_ymin_label = new System.Windows.Forms.Label();
            this.pot_ymin_val = new System.Windows.Forms.TextBox();
            this.pot_xmax_label = new System.Windows.Forms.Label();
            this.pot_xmax_val = new System.Windows.Forms.TextBox();
            this.pot_plot_label = new System.Windows.Forms.Label();
            this.pot_xmin_label = new System.Windows.Forms.Label();
            this.pot_xmin_val = new System.Windows.Forms.TextBox();
            this.dens_ymax_label = new System.Windows.Forms.Label();
            this.dens_ymax_val = new System.Windows.Forms.TextBox();
            this.dens_ymin_label = new System.Windows.Forms.Label();
            this.dens_ymin_val = new System.Windows.Forms.TextBox();
            this.dens_xmax_label = new System.Windows.Forms.Label();
            this.dens_xmax_val = new System.Windows.Forms.TextBox();
            this.dens_plot_label = new System.Windows.Forms.Label();
            this.dens_xmin_label = new System.Windows.Forms.Label();
            this.dens_xmin_val = new System.Windows.Forms.TextBox();
            this.run_button = new System.Windows.Forms.Button();
            this.temperature_val_label = new System.Windows.Forms.Label();
            this.temperature_label = new System.Windows.Forms.Label();
            this.count_no_label = new System.Windows.Forms.Label();
            this.refresh_button = new System.Windows.Forms.Button();
            this.count_label = new System.Windows.Forms.Label();
            this.step_button = new System.Windows.Forms.Button();
            this.group1Dboundaries = new System.Windows.Forms.GroupBox();
            this.boundary1Ddescriptor_label = new System.Windows.Forms.Label();
            this.bottomV1Dval_label = new System.Windows.Forms.Label();
            this.topV1Dval = new System.Windows.Forms.TextBox();
            this.topV1Dval_label = new System.Windows.Forms.Label();
            this.bottomV1Dval = new System.Windows.Forms.TextBox();
            this.groupConversion = new System.Windows.Forms.GroupBox();
            this.ny1Dval_label = new System.Windows.Forms.Label();
            this.nx1Dval_label = new System.Windows.Forms.Label();
            this.ny1Dval = new System.Windows.Forms.TextBox();
            this.nx1Dval = new System.Windows.Forms.TextBox();
            this.groupPhysical = new System.Windows.Forms.GroupBox();
            this.temperature1Dval_label = new System.Windows.Forms.Label();
            this.temperature1Dval = new System.Windows.Forms.TextBox();
            this.groupDFT = new System.Windows.Forms.GroupBox();
            this.dftzmin1Dval_label = new System.Windows.Forms.Label();
            this.dftzmin1Dval = new System.Windows.Forms.TextBox();
            this.dftnz1Dval_label = new System.Windows.Forms.Label();
            this.dft1DCheck = new System.Windows.Forms.CheckBox();
            this.dftnz1Dval = new System.Windows.Forms.TextBox();
            this.groupPotential = new System.Windows.Forms.GroupBox();
            this.nz1Dval_label = new System.Windows.Forms.Label();
            this.nz1Dval = new System.Windows.Forms.TextBox();
            this.dz1Dval_label = new System.Windows.Forms.Label();
            this.dz1Dval = new System.Windows.Forms.TextBox();
            this.dopentDescriptor = new System.Windows.Forms.Label();
            this.otherinputs = new System.Windows.Forms.TabPage();
            this.groupConvergence = new System.Windows.Forms.GroupBox();
            this.max_iterations1d_val = new System.Windows.Forms.TextBox();
            this.max_iterations1d_labal = new System.Windows.Forms.Label();
            this.dopent_convergenceparams_label = new System.Windows.Forms.Label();
            this.main_convergenceparams_label = new System.Windows.Forms.Label();
            this.maxiteration_val = new System.Windows.Forms.TextBox();
            this.max_iteration_label = new System.Windows.Forms.Label();
            this.toleranceval = new System.Windows.Forms.TextBox();
            this.tolerance_label = new System.Windows.Forms.Label();
            this.tolerance1Dval = new System.Windows.Forms.TextBox();
            this.tolerance1D_label = new System.Windows.Forms.Label();
            this.dimensionLabel = new System.Windows.Forms.Label();
            this.dimensionality = new System.Windows.Forms.ListBox();
            this.includebatchdata = new System.Windows.Forms.CheckBox();
            this.checkpointing = new System.Windows.Forms.CheckBox();
            this.subsplit = new System.Windows.Forms.SplitContainer();
            this.conduction_band = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.density = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.carrier_dopent_density_Text = new System.Windows.Forms.TextBox();
            this.carrier_density_dopent_label = new System.Windows.Forms.Label();
            this.carrier_density_dopent_units = new System.Windows.Forms.Label();
            this.menuStrip1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer)).BeginInit();
            this.splitContainer.Panel1.SuspendLayout();
            this.splitContainer.Panel2.SuspendLayout();
            this.splitContainer.SuspendLayout();
            this.inputTab.SuspendLayout();
            this.bandstructure.SuspendLayout();
            this.group_newlayer.SuspendLayout();
            this.densityinputs.SuspendLayout();
            this.dopentinputs.SuspendLayout();
            this.groupPlotting.SuspendLayout();
            this.group1Dboundaries.SuspendLayout();
            this.groupConversion.SuspendLayout();
            this.groupPhysical.SuspendLayout();
            this.groupDFT.SuspendLayout();
            this.groupPotential.SuspendLayout();
            this.otherinputs.SuspendLayout();
            this.groupConvergence.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.subsplit)).BeginInit();
            this.subsplit.Panel1.SuspendLayout();
            this.subsplit.Panel2.SuspendLayout();
            this.subsplit.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.conduction_band)).BeginInit();
            ((System.ComponentModel.ISupportInitialize)(this.density)).BeginInit();
            this.SuspendLayout();
            // 
            // menuStrip1
            // 
            this.menuStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.fileToolStripMenuItem,
            this.viewToolStripMenuItem,
            this.settingToolStripMenuItem,
            this.helpToolStripMenuItem});
            this.menuStrip1.Location = new System.Drawing.Point(0, 0);
            this.menuStrip1.Name = "menuStrip1";
            this.menuStrip1.Size = new System.Drawing.Size(1416, 24);
            this.menuStrip1.TabIndex = 0;
            this.menuStrip1.Text = "menuStrip1";
            // 
            // fileToolStripMenuItem
            // 
            this.fileToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.mnuNew,
            this.mnuOpen,
            this.toolStripMenuItem1,
            this.mnuSave,
            this.mnuSaveAs,
            this.toolStripMenuItem2,
            this.mnuExit});
            this.fileToolStripMenuItem.Name = "fileToolStripMenuItem";
            this.fileToolStripMenuItem.Size = new System.Drawing.Size(37, 20);
            this.fileToolStripMenuItem.Text = "File";
            // 
            // mnuNew
            // 
            this.mnuNew.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.mnuNewBandStructure,
            this.mnuNewInputFile,
            this.solverConfigToolStripMenuItem});
            this.mnuNew.Name = "mnuNew";
            this.mnuNew.Size = new System.Drawing.Size(138, 22);
            this.mnuNew.Text = "New";
            // 
            // mnuNewBandStructure
            // 
            this.mnuNewBandStructure.Name = "mnuNewBandStructure";
            this.mnuNewBandStructure.Size = new System.Drawing.Size(152, 22);
            this.mnuNewBandStructure.Text = "Band Structure";
            // 
            // mnuNewInputFile
            // 
            this.mnuNewInputFile.Name = "mnuNewInputFile";
            this.mnuNewInputFile.Size = new System.Drawing.Size(152, 22);
            this.mnuNewInputFile.Text = "Input File";
            // 
            // solverConfigToolStripMenuItem
            // 
            this.solverConfigToolStripMenuItem.Name = "solverConfigToolStripMenuItem";
            this.solverConfigToolStripMenuItem.Size = new System.Drawing.Size(152, 22);
            this.solverConfigToolStripMenuItem.Text = "Solver Config";
            // 
            // mnuOpen
            // 
            this.mnuOpen.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.mnuOpenBandStructure,
            this.mnuOpenInputFile,
            this.solverConfigToolStripMenuItem1});
            this.mnuOpen.Name = "mnuOpen";
            this.mnuOpen.Size = new System.Drawing.Size(138, 22);
            this.mnuOpen.Text = "Open";
            // 
            // mnuOpenBandStructure
            // 
            this.mnuOpenBandStructure.Name = "mnuOpenBandStructure";
            this.mnuOpenBandStructure.Size = new System.Drawing.Size(152, 22);
            this.mnuOpenBandStructure.Text = "Band Structure";
            this.mnuOpenBandStructure.Click += new System.EventHandler(this.mnuOpenBandStructure_Click);
            // 
            // mnuOpenInputFile
            // 
            this.mnuOpenInputFile.Name = "mnuOpenInputFile";
            this.mnuOpenInputFile.Size = new System.Drawing.Size(152, 22);
            this.mnuOpenInputFile.Text = "Input File";
            // 
            // solverConfigToolStripMenuItem1
            // 
            this.solverConfigToolStripMenuItem1.Name = "solverConfigToolStripMenuItem1";
            this.solverConfigToolStripMenuItem1.Size = new System.Drawing.Size(152, 22);
            this.solverConfigToolStripMenuItem1.Text = "Solver Config";
            // 
            // toolStripMenuItem1
            // 
            this.toolStripMenuItem1.Name = "toolStripMenuItem1";
            this.toolStripMenuItem1.Size = new System.Drawing.Size(135, 6);
            // 
            // mnuSave
            // 
            this.mnuSave.Name = "mnuSave";
            this.mnuSave.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.S)));
            this.mnuSave.Size = new System.Drawing.Size(138, 22);
            this.mnuSave.Text = "Save";
            // 
            // mnuSaveAs
            // 
            this.mnuSaveAs.Name = "mnuSaveAs";
            this.mnuSaveAs.Size = new System.Drawing.Size(138, 22);
            this.mnuSaveAs.Text = "Save As";
            // 
            // toolStripMenuItem2
            // 
            this.toolStripMenuItem2.Name = "toolStripMenuItem2";
            this.toolStripMenuItem2.Size = new System.Drawing.Size(135, 6);
            // 
            // mnuExit
            // 
            this.mnuExit.Name = "mnuExit";
            this.mnuExit.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.Q)));
            this.mnuExit.Size = new System.Drawing.Size(138, 22);
            this.mnuExit.Text = "Exit";
            this.mnuExit.Click += new System.EventHandler(this.exitToolStripMenuItem_Click);
            // 
            // viewToolStripMenuItem
            // 
            this.viewToolStripMenuItem.Name = "viewToolStripMenuItem";
            this.viewToolStripMenuItem.Size = new System.Drawing.Size(44, 20);
            this.viewToolStripMenuItem.Text = "View";
            // 
            // settingToolStripMenuItem
            // 
            this.settingToolStripMenuItem.Name = "settingToolStripMenuItem";
            this.settingToolStripMenuItem.Size = new System.Drawing.Size(61, 20);
            this.settingToolStripMenuItem.Text = "Options";
            // 
            // helpToolStripMenuItem
            // 
            this.helpToolStripMenuItem.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.gettingStartedToolStripMenuItem,
            this.aboutToolStripMenuItem});
            this.helpToolStripMenuItem.Name = "helpToolStripMenuItem";
            this.helpToolStripMenuItem.Size = new System.Drawing.Size(44, 20);
            this.helpToolStripMenuItem.Text = "Help";
            // 
            // gettingStartedToolStripMenuItem
            // 
            this.gettingStartedToolStripMenuItem.Name = "gettingStartedToolStripMenuItem";
            this.gettingStartedToolStripMenuItem.Size = new System.Drawing.Size(153, 22);
            this.gettingStartedToolStripMenuItem.Text = "Getting Started";
            // 
            // aboutToolStripMenuItem
            // 
            this.aboutToolStripMenuItem.Name = "aboutToolStripMenuItem";
            this.aboutToolStripMenuItem.Size = new System.Drawing.Size(153, 22);
            this.aboutToolStripMenuItem.Text = "About";
            // 
            // splitContainer
            // 
            this.splitContainer.Dock = System.Windows.Forms.DockStyle.Fill;
            this.splitContainer.Location = new System.Drawing.Point(0, 24);
            this.splitContainer.Name = "splitContainer";
            // 
            // splitContainer.Panel1
            // 
            this.splitContainer.Panel1.AccessibleRole = System.Windows.Forms.AccessibleRole.OutlineButton;
            this.splitContainer.Panel1.Controls.Add(this.inputTab);
            // 
            // splitContainer.Panel2
            // 
            this.splitContainer.Panel2.Controls.Add(this.subsplit);
            this.splitContainer.Size = new System.Drawing.Size(1416, 649);
            this.splitContainer.SplitterDistance = 734;
            this.splitContainer.TabIndex = 1;
            // 
            // inputTab
            // 
            this.inputTab.Controls.Add(this.bandstructure);
            this.inputTab.Controls.Add(this.potentialinputs);
            this.inputTab.Controls.Add(this.densityinputs);
            this.inputTab.Controls.Add(this.dopentinputs);
            this.inputTab.Controls.Add(this.otherinputs);
            this.inputTab.Location = new System.Drawing.Point(3, 3);
            this.inputTab.Name = "inputTab";
            this.inputTab.SelectedIndex = 0;
            this.inputTab.Size = new System.Drawing.Size(729, 646);
            this.inputTab.TabIndex = 0;
            // 
            // bandstructure
            // 
            this.bandstructure.Controls.Add(this.total_thickness_label);
            this.bandstructure.Controls.Add(this.total_thickness_val);
            this.bandstructure.Controls.Add(this.bandstructurefilename);
            this.bandstructure.Controls.Add(this.bandstructurefile_label);
            this.bandstructure.Controls.Add(this.bandstructure_create_button);
            this.bandstructure.Controls.Add(this.group_newlayer);
            this.bandstructure.Controls.Add(this.layerselect_label);
            this.bandstructure.Controls.Add(this.bandstructure_list);
            this.bandstructure.Controls.Add(this.editLayer_Button);
            this.bandstructure.Controls.Add(this.deleteLayer_Button);
            this.bandstructure.Controls.Add(this.addlayer_Button);
            this.bandstructure.Controls.Add(this.bandstructureCombo_label);
            this.bandstructure.Controls.Add(this.bandstructureCombo);
            this.bandstructure.Location = new System.Drawing.Point(4, 22);
            this.bandstructure.Name = "bandstructure";
            this.bandstructure.Padding = new System.Windows.Forms.Padding(3);
            this.bandstructure.Size = new System.Drawing.Size(721, 620);
            this.bandstructure.TabIndex = 0;
            this.bandstructure.Text = "Band Structure";
            this.bandstructure.UseVisualStyleBackColor = true;
            // 
            // total_thickness_label
            // 
            this.total_thickness_label.AutoSize = true;
            this.total_thickness_label.Location = new System.Drawing.Point(15, 526);
            this.total_thickness_label.Name = "total_thickness_label";
            this.total_thickness_label.Size = new System.Drawing.Size(83, 13);
            this.total_thickness_label.TabIndex = 12;
            this.total_thickness_label.Text = "Total Thickness";
            // 
            // total_thickness_val
            // 
            this.total_thickness_val.Location = new System.Drawing.Point(100, 523);
            this.total_thickness_val.Name = "total_thickness_val";
            this.total_thickness_val.ReadOnly = true;
            this.total_thickness_val.Size = new System.Drawing.Size(92, 20);
            this.total_thickness_val.TabIndex = 11;
            this.total_thickness_val.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // bandstructurefilename
            // 
            this.bandstructurefilename.Location = new System.Drawing.Point(32, 584);
            this.bandstructurefilename.Name = "bandstructurefilename";
            this.bandstructurefilename.Size = new System.Drawing.Size(294, 20);
            this.bandstructurefilename.TabIndex = 10;
            this.bandstructurefilename.Text = "BandStructure.txt";
            // 
            // bandstructurefile_label
            // 
            this.bandstructurefile_label.AutoSize = true;
            this.bandstructurefile_label.Location = new System.Drawing.Point(29, 568);
            this.bandstructurefile_label.Name = "bandstructurefile_label";
            this.bandstructurefile_label.Size = new System.Drawing.Size(128, 13);
            this.bandstructurefile_label.TabIndex = 9;
            this.bandstructurefile_label.Text = "Band Structure File Name";
            // 
            // bandstructure_create_button
            // 
            this.bandstructure_create_button.Location = new System.Drawing.Point(342, 578);
            this.bandstructure_create_button.Name = "bandstructure_create_button";
            this.bandstructure_create_button.Size = new System.Drawing.Size(143, 31);
            this.bandstructure_create_button.TabIndex = 8;
            this.bandstructure_create_button.Text = "Create Band Structure File";
            this.bandstructure_create_button.UseVisualStyleBackColor = true;
            this.bandstructure_create_button.Click += new System.EventHandler(this.bandstructure_create_button_Click);
            // 
            // group_newlayer
            // 
            this.group_newlayer.Controls.Add(this.newlayer_naval);
            this.group_newlayer.Controls.Add(this.newlayer_na_label);
            this.group_newlayer.Controls.Add(this.newlayer_ndval);
            this.group_newlayer.Controls.Add(this.newlayer_nd_label);
            this.group_newlayer.Controls.Add(this.newlayer_xval);
            this.group_newlayer.Controls.Add(this.newlayer_x_label);
            this.group_newlayer.Controls.Add(this.newlayer_thicknessval);
            this.group_newlayer.Controls.Add(this.newlayer_thickness_label);
            this.group_newlayer.Controls.Add(this.material_label);
            this.group_newlayer.Controls.Add(this.material_combo);
            this.group_newlayer.Location = new System.Drawing.Point(14, 95);
            this.group_newlayer.Name = "group_newlayer";
            this.group_newlayer.Size = new System.Drawing.Size(574, 71);
            this.group_newlayer.TabIndex = 7;
            this.group_newlayer.TabStop = false;
            this.group_newlayer.Text = "New Layer";
            // 
            // newlayer_naval
            // 
            this.newlayer_naval.Location = new System.Drawing.Point(463, 40);
            this.newlayer_naval.Name = "newlayer_naval";
            this.newlayer_naval.Size = new System.Drawing.Size(100, 20);
            this.newlayer_naval.TabIndex = 9;
            // 
            // newlayer_na_label
            // 
            this.newlayer_na_label.AutoSize = true;
            this.newlayer_na_label.Location = new System.Drawing.Point(462, 24);
            this.newlayer_na_label.Name = "newlayer_na_label";
            this.newlayer_na_label.Size = new System.Drawing.Size(61, 13);
            this.newlayer_na_label.TabIndex = 8;
            this.newlayer_na_label.Text = "Na / cm^-3";
            // 
            // newlayer_ndval
            // 
            this.newlayer_ndval.Location = new System.Drawing.Point(357, 40);
            this.newlayer_ndval.Name = "newlayer_ndval";
            this.newlayer_ndval.Size = new System.Drawing.Size(100, 20);
            this.newlayer_ndval.TabIndex = 7;
            // 
            // newlayer_nd_label
            // 
            this.newlayer_nd_label.AutoSize = true;
            this.newlayer_nd_label.Location = new System.Drawing.Point(356, 24);
            this.newlayer_nd_label.Name = "newlayer_nd_label";
            this.newlayer_nd_label.Size = new System.Drawing.Size(61, 13);
            this.newlayer_nd_label.TabIndex = 6;
            this.newlayer_nd_label.Text = "Nd / cm^-3";
            // 
            // newlayer_xval
            // 
            this.newlayer_xval.Enabled = false;
            this.newlayer_xval.Location = new System.Drawing.Point(251, 40);
            this.newlayer_xval.Name = "newlayer_xval";
            this.newlayer_xval.Size = new System.Drawing.Size(100, 20);
            this.newlayer_xval.TabIndex = 5;
            // 
            // newlayer_x_label
            // 
            this.newlayer_x_label.AutoSize = true;
            this.newlayer_x_label.Location = new System.Drawing.Point(250, 24);
            this.newlayer_x_label.Name = "newlayer_x_label";
            this.newlayer_x_label.Size = new System.Drawing.Size(103, 13);
            this.newlayer_x_label.TabIndex = 4;
            this.newlayer_x_label.Text = "Alloy Composition (x)";
            // 
            // newlayer_thicknessval
            // 
            this.newlayer_thicknessval.Location = new System.Drawing.Point(145, 40);
            this.newlayer_thicknessval.Name = "newlayer_thicknessval";
            this.newlayer_thicknessval.Size = new System.Drawing.Size(100, 20);
            this.newlayer_thicknessval.TabIndex = 3;
            // 
            // newlayer_thickness_label
            // 
            this.newlayer_thickness_label.AutoSize = true;
            this.newlayer_thickness_label.Location = new System.Drawing.Point(144, 24);
            this.newlayer_thickness_label.Name = "newlayer_thickness_label";
            this.newlayer_thickness_label.Size = new System.Drawing.Size(81, 13);
            this.newlayer_thickness_label.TabIndex = 2;
            this.newlayer_thickness_label.Text = "Thickness / nm";
            // 
            // material_label
            // 
            this.material_label.AutoSize = true;
            this.material_label.Location = new System.Drawing.Point(15, 24);
            this.material_label.Name = "material_label";
            this.material_label.Size = new System.Drawing.Size(44, 13);
            this.material_label.TabIndex = 1;
            this.material_label.Text = "Material";
            // 
            // material_combo
            // 
            this.material_combo.FormattingEnabled = true;
            this.material_combo.Items.AddRange(new object[] {
            "GaAs",
            "AlGaAs",
            "PMMA",
            "InAlAs",
            "InGaAs",
            "Metal",
            "Air",
            "Composite"});
            this.material_combo.Location = new System.Drawing.Point(18, 40);
            this.material_combo.Name = "material_combo";
            this.material_combo.Size = new System.Drawing.Size(121, 21);
            this.material_combo.TabIndex = 0;
            this.material_combo.Text = "GaAs";
            this.material_combo.SelectedIndexChanged += new System.EventHandler(this.material_combo_SelectedIndexChanged);
            // 
            // layerselect_label
            // 
            this.layerselect_label.AutoSize = true;
            this.layerselect_label.Location = new System.Drawing.Point(18, 22);
            this.layerselect_label.Name = "layerselect_label";
            this.layerselect_label.Size = new System.Drawing.Size(80, 13);
            this.layerselect_label.TabIndex = 6;
            this.layerselect_label.Text = "Layer Selection";
            // 
            // bandstructure_list
            // 
            this.bandstructure_list.Columns.AddRange(new System.Windows.Forms.ColumnHeader[] {
            this.material_header,
            this.thickness_header,
            this.alloy_header,
            this.donor_header,
            this.acceptor_header});
            this.bandstructure_list.ForeColor = System.Drawing.SystemColors.WindowText;
            this.bandstructure_list.GridLines = true;
            listViewItem1.StateImageIndex = 0;
            listViewItem1.UseItemStyleForSubItems = false;
            listViewItem2.UseItemStyleForSubItems = false;
            this.bandstructure_list.Items.AddRange(new System.Windows.Forms.ListViewItem[] {
            listViewItem1,
            listViewItem2});
            this.bandstructure_list.Location = new System.Drawing.Point(32, 198);
            this.bandstructure_list.Name = "bandstructure_list";
            this.bandstructure_list.Size = new System.Drawing.Size(448, 319);
            this.bandstructure_list.TabIndex = 5;
            this.bandstructure_list.UseCompatibleStateImageBehavior = false;
            this.bandstructure_list.View = System.Windows.Forms.View.Details;
            // 
            // material_header
            // 
            this.material_header.Text = "Material";
            this.material_header.Width = 65;
            // 
            // thickness_header
            // 
            this.thickness_header.Text = "Thickness / nm";
            this.thickness_header.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.thickness_header.Width = 90;
            // 
            // alloy_header
            // 
            this.alloy_header.Text = "Alloy Composition (x)";
            this.alloy_header.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.alloy_header.Width = 112;
            // 
            // donor_header
            // 
            this.donor_header.Text = "Nd / cm^-3";
            this.donor_header.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.donor_header.Width = 74;
            // 
            // acceptor_header
            // 
            this.acceptor_header.Text = "Na / cm^-3";
            this.acceptor_header.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            this.acceptor_header.Width = 103;
            // 
            // editLayer_Button
            // 
            this.editLayer_Button.Location = new System.Drawing.Point(342, 37);
            this.editLayer_Button.Name = "editLayer_Button";
            this.editLayer_Button.Size = new System.Drawing.Size(91, 22);
            this.editLayer_Button.TabIndex = 4;
            this.editLayer_Button.Text = "Edit Layer";
            this.editLayer_Button.UseVisualStyleBackColor = true;
            this.editLayer_Button.Click += new System.EventHandler(this.editLayer_Button_Click);
            // 
            // deleteLayer_Button
            // 
            this.deleteLayer_Button.Location = new System.Drawing.Point(245, 37);
            this.deleteLayer_Button.Name = "deleteLayer_Button";
            this.deleteLayer_Button.Size = new System.Drawing.Size(91, 22);
            this.deleteLayer_Button.TabIndex = 3;
            this.deleteLayer_Button.Text = "Delete Layer";
            this.deleteLayer_Button.UseVisualStyleBackColor = true;
            this.deleteLayer_Button.Click += new System.EventHandler(this.deleteLayer_Button_Click);
            // 
            // addlayer_Button
            // 
            this.addlayer_Button.Location = new System.Drawing.Point(148, 37);
            this.addlayer_Button.Name = "addlayer_Button";
            this.addlayer_Button.Size = new System.Drawing.Size(91, 22);
            this.addlayer_Button.TabIndex = 2;
            this.addlayer_Button.Text = "Add Layer";
            this.addlayer_Button.UseVisualStyleBackColor = true;
            this.addlayer_Button.Click += new System.EventHandler(this.addlayer_Button_Click);
            // 
            // bandstructureCombo_label
            // 
            this.bandstructureCombo_label.AutoSize = true;
            this.bandstructureCombo_label.Location = new System.Drawing.Point(29, 182);
            this.bandstructureCombo_label.Name = "bandstructureCombo_label";
            this.bandstructureCombo_label.Size = new System.Drawing.Size(115, 13);
            this.bandstructureCombo_label.TabIndex = 1;
            this.bandstructureCombo_label.Text = "Current Band Structure";
            // 
            // bandstructureCombo
            // 
            this.bandstructureCombo.FormattingEnabled = true;
            this.bandstructureCombo.Items.AddRange(new object[] {
            "Surface",
            "Substrate"});
            this.bandstructureCombo.Location = new System.Drawing.Point(21, 38);
            this.bandstructureCombo.Name = "bandstructureCombo";
            this.bandstructureCombo.Size = new System.Drawing.Size(121, 21);
            this.bandstructureCombo.TabIndex = 0;
            // 
            // potentialinputs
            // 
            this.potentialinputs.Location = new System.Drawing.Point(4, 22);
            this.potentialinputs.Name = "potentialinputs";
            this.potentialinputs.Padding = new System.Windows.Forms.Padding(3);
            this.potentialinputs.Size = new System.Drawing.Size(721, 620);
            this.potentialinputs.TabIndex = 1;
            this.potentialinputs.Text = "Potential";
            this.potentialinputs.UseVisualStyleBackColor = true;
            // 
            // densityinputs
            // 
            this.densityinputs.Controls.Add(this.output_label);
            this.densityinputs.Location = new System.Drawing.Point(4, 22);
            this.densityinputs.Name = "densityinputs";
            this.densityinputs.Size = new System.Drawing.Size(721, 620);
            this.densityinputs.TabIndex = 2;
            this.densityinputs.Text = "Density";
            this.densityinputs.UseVisualStyleBackColor = true;
            // 
            // output_label
            // 
            this.output_label.AutoSize = true;
            this.output_label.Location = new System.Drawing.Point(38, 410);
            this.output_label.Name = "output_label";
            this.output_label.Size = new System.Drawing.Size(0, 13);
            this.output_label.TabIndex = 3;
            // 
            // dopentinputs
            // 
            this.dopentinputs.Controls.Add(this.carrier_density_dopent_units);
            this.dopentinputs.Controls.Add(this.carrier_density_dopent_label);
            this.dopentinputs.Controls.Add(this.carrier_dopent_density_Text);
            this.dopentinputs.Controls.Add(this.groupPlotting);
            this.dopentinputs.Controls.Add(this.run_button);
            this.dopentinputs.Controls.Add(this.temperature_val_label);
            this.dopentinputs.Controls.Add(this.temperature_label);
            this.dopentinputs.Controls.Add(this.count_no_label);
            this.dopentinputs.Controls.Add(this.refresh_button);
            this.dopentinputs.Controls.Add(this.count_label);
            this.dopentinputs.Controls.Add(this.step_button);
            this.dopentinputs.Controls.Add(this.group1Dboundaries);
            this.dopentinputs.Controls.Add(this.groupConversion);
            this.dopentinputs.Controls.Add(this.groupPhysical);
            this.dopentinputs.Controls.Add(this.groupDFT);
            this.dopentinputs.Controls.Add(this.groupPotential);
            this.dopentinputs.Controls.Add(this.dopentDescriptor);
            this.dopentinputs.Location = new System.Drawing.Point(4, 22);
            this.dopentinputs.Name = "dopentinputs";
            this.dopentinputs.Size = new System.Drawing.Size(721, 620);
            this.dopentinputs.TabIndex = 3;
            this.dopentinputs.Text = "Dopents";
            this.dopentinputs.UseVisualStyleBackColor = true;
            // 
            // groupPlotting
            // 
            this.groupPlotting.Controls.Add(this.pot_ymax_label);
            this.groupPlotting.Controls.Add(this.pot_ymax_val);
            this.groupPlotting.Controls.Add(this.pot_ymin_label);
            this.groupPlotting.Controls.Add(this.pot_ymin_val);
            this.groupPlotting.Controls.Add(this.pot_xmax_label);
            this.groupPlotting.Controls.Add(this.pot_xmax_val);
            this.groupPlotting.Controls.Add(this.pot_plot_label);
            this.groupPlotting.Controls.Add(this.pot_xmin_label);
            this.groupPlotting.Controls.Add(this.pot_xmin_val);
            this.groupPlotting.Controls.Add(this.dens_ymax_label);
            this.groupPlotting.Controls.Add(this.dens_ymax_val);
            this.groupPlotting.Controls.Add(this.dens_ymin_label);
            this.groupPlotting.Controls.Add(this.dens_ymin_val);
            this.groupPlotting.Controls.Add(this.dens_xmax_label);
            this.groupPlotting.Controls.Add(this.dens_xmax_val);
            this.groupPlotting.Controls.Add(this.dens_plot_label);
            this.groupPlotting.Controls.Add(this.dens_xmin_label);
            this.groupPlotting.Controls.Add(this.dens_xmin_val);
            this.groupPlotting.Location = new System.Drawing.Point(304, 205);
            this.groupPlotting.Name = "groupPlotting";
            this.groupPlotting.Size = new System.Drawing.Size(392, 154);
            this.groupPlotting.TabIndex = 13;
            this.groupPlotting.TabStop = false;
            this.groupPlotting.Text = "Plotting Parameters";
            // 
            // pot_ymax_label
            // 
            this.pot_ymax_label.AutoSize = true;
            this.pot_ymax_label.Location = new System.Drawing.Point(306, 120);
            this.pot_ymax_label.Name = "pot_ymax_label";
            this.pot_ymax_label.Size = new System.Drawing.Size(77, 13);
            this.pot_ymax_label.TabIndex = 29;
            this.pot_ymax_label.Text = "ymax / nm^-1?";
            // 
            // pot_ymax_val
            // 
            this.pot_ymax_val.Location = new System.Drawing.Point(206, 117);
            this.pot_ymax_val.Name = "pot_ymax_val";
            this.pot_ymax_val.Size = new System.Drawing.Size(94, 20);
            this.pot_ymax_val.TabIndex = 28;
            // 
            // pot_ymin_label
            // 
            this.pot_ymin_label.AutoSize = true;
            this.pot_ymin_label.Location = new System.Drawing.Point(306, 94);
            this.pot_ymin_label.Name = "pot_ymin_label";
            this.pot_ymin_label.Size = new System.Drawing.Size(74, 13);
            this.pot_ymin_label.TabIndex = 27;
            this.pot_ymin_label.Text = "ymin / nm^-1?";
            // 
            // pot_ymin_val
            // 
            this.pot_ymin_val.Location = new System.Drawing.Point(206, 91);
            this.pot_ymin_val.Name = "pot_ymin_val";
            this.pot_ymin_val.Size = new System.Drawing.Size(94, 20);
            this.pot_ymin_val.TabIndex = 26;
            // 
            // pot_xmax_label
            // 
            this.pot_xmax_label.AutoSize = true;
            this.pot_xmax_label.Location = new System.Drawing.Point(306, 68);
            this.pot_xmax_label.Name = "pot_xmax_label";
            this.pot_xmax_label.Size = new System.Drawing.Size(77, 13);
            this.pot_xmax_label.TabIndex = 25;
            this.pot_xmax_label.Text = "xmax / nm^-1?";
            // 
            // pot_xmax_val
            // 
            this.pot_xmax_val.Location = new System.Drawing.Point(206, 65);
            this.pot_xmax_val.Name = "pot_xmax_val";
            this.pot_xmax_val.Size = new System.Drawing.Size(94, 20);
            this.pot_xmax_val.TabIndex = 24;
            // 
            // pot_plot_label
            // 
            this.pot_plot_label.AutoSize = true;
            this.pot_plot_label.Location = new System.Drawing.Point(203, 23);
            this.pot_plot_label.Name = "pot_plot_label";
            this.pot_plot_label.Size = new System.Drawing.Size(48, 13);
            this.pot_plot_label.TabIndex = 23;
            this.pot_plot_label.Text = "Potential";
            // 
            // pot_xmin_label
            // 
            this.pot_xmin_label.AutoSize = true;
            this.pot_xmin_label.Location = new System.Drawing.Point(306, 42);
            this.pot_xmin_label.Name = "pot_xmin_label";
            this.pot_xmin_label.Size = new System.Drawing.Size(74, 13);
            this.pot_xmin_label.TabIndex = 22;
            this.pot_xmin_label.Text = "xmin / nm^-1?";
            // 
            // pot_xmin_val
            // 
            this.pot_xmin_val.Location = new System.Drawing.Point(206, 39);
            this.pot_xmin_val.Name = "pot_xmin_val";
            this.pot_xmin_val.Size = new System.Drawing.Size(94, 20);
            this.pot_xmin_val.TabIndex = 21;
            // 
            // dens_ymax_label
            // 
            this.dens_ymax_label.AutoSize = true;
            this.dens_ymax_label.Location = new System.Drawing.Point(117, 120);
            this.dens_ymax_label.Name = "dens_ymax_label";
            this.dens_ymax_label.Size = new System.Drawing.Size(77, 13);
            this.dens_ymax_label.TabIndex = 20;
            this.dens_ymax_label.Text = "ymax / nm^-1?";
            // 
            // dens_ymax_val
            // 
            this.dens_ymax_val.Location = new System.Drawing.Point(17, 117);
            this.dens_ymax_val.Name = "dens_ymax_val";
            this.dens_ymax_val.Size = new System.Drawing.Size(94, 20);
            this.dens_ymax_val.TabIndex = 19;
            // 
            // dens_ymin_label
            // 
            this.dens_ymin_label.AutoSize = true;
            this.dens_ymin_label.Location = new System.Drawing.Point(117, 94);
            this.dens_ymin_label.Name = "dens_ymin_label";
            this.dens_ymin_label.Size = new System.Drawing.Size(74, 13);
            this.dens_ymin_label.TabIndex = 18;
            this.dens_ymin_label.Text = "ymin / nm^-1?";
            // 
            // dens_ymin_val
            // 
            this.dens_ymin_val.Location = new System.Drawing.Point(17, 91);
            this.dens_ymin_val.Name = "dens_ymin_val";
            this.dens_ymin_val.Size = new System.Drawing.Size(94, 20);
            this.dens_ymin_val.TabIndex = 17;
            // 
            // dens_xmax_label
            // 
            this.dens_xmax_label.AutoSize = true;
            this.dens_xmax_label.Location = new System.Drawing.Point(117, 68);
            this.dens_xmax_label.Name = "dens_xmax_label";
            this.dens_xmax_label.Size = new System.Drawing.Size(77, 13);
            this.dens_xmax_label.TabIndex = 16;
            this.dens_xmax_label.Text = "xmax / nm^-1?";
            // 
            // dens_xmax_val
            // 
            this.dens_xmax_val.Location = new System.Drawing.Point(17, 65);
            this.dens_xmax_val.Name = "dens_xmax_val";
            this.dens_xmax_val.Size = new System.Drawing.Size(94, 20);
            this.dens_xmax_val.TabIndex = 15;
            // 
            // dens_plot_label
            // 
            this.dens_plot_label.AutoSize = true;
            this.dens_plot_label.Location = new System.Drawing.Point(14, 23);
            this.dens_plot_label.Name = "dens_plot_label";
            this.dens_plot_label.Size = new System.Drawing.Size(42, 13);
            this.dens_plot_label.TabIndex = 14;
            this.dens_plot_label.Text = "Density";
            // 
            // dens_xmin_label
            // 
            this.dens_xmin_label.AutoSize = true;
            this.dens_xmin_label.Location = new System.Drawing.Point(117, 42);
            this.dens_xmin_label.Name = "dens_xmin_label";
            this.dens_xmin_label.Size = new System.Drawing.Size(74, 13);
            this.dens_xmin_label.TabIndex = 5;
            this.dens_xmin_label.Text = "xmin / nm^-1?";
            // 
            // dens_xmin_val
            // 
            this.dens_xmin_val.Location = new System.Drawing.Point(17, 39);
            this.dens_xmin_val.Name = "dens_xmin_val";
            this.dens_xmin_val.Size = new System.Drawing.Size(94, 20);
            this.dens_xmin_val.TabIndex = 4;
            // 
            // run_button
            // 
            this.run_button.Location = new System.Drawing.Point(508, 60);
            this.run_button.Name = "run_button";
            this.run_button.Size = new System.Drawing.Size(97, 82);
            this.run_button.TabIndex = 12;
            this.run_button.Text = "Run";
            this.run_button.UseVisualStyleBackColor = true;
            this.run_button.Click += new System.EventHandler(this.run_button_Click);
            // 
            // temperature_val_label
            // 
            this.temperature_val_label.AutoSize = true;
            this.temperature_val_label.Location = new System.Drawing.Point(429, 173);
            this.temperature_val_label.Name = "temperature_val_label";
            this.temperature_val_label.Size = new System.Drawing.Size(27, 13);
            this.temperature_val_label.TabIndex = 11;
            this.temperature_val_label.Text = "N/A";
            // 
            // temperature_label
            // 
            this.temperature_label.AutoSize = true;
            this.temperature_label.Location = new System.Drawing.Point(353, 173);
            this.temperature_label.Name = "temperature_label";
            this.temperature_label.Size = new System.Drawing.Size(79, 13);
            this.temperature_label.TabIndex = 10;
            this.temperature_label.Text = "Temperature = ";
            // 
            // count_no_label
            // 
            this.count_no_label.AutoSize = true;
            this.count_no_label.Location = new System.Drawing.Point(396, 156);
            this.count_no_label.Name = "count_no_label";
            this.count_no_label.Size = new System.Drawing.Size(13, 13);
            this.count_no_label.TabIndex = 9;
            this.count_no_label.Text = "0";
            // 
            // refresh_button
            // 
            this.refresh_button.Location = new System.Drawing.Point(508, 159);
            this.refresh_button.Name = "refresh_button";
            this.refresh_button.Size = new System.Drawing.Size(75, 23);
            this.refresh_button.TabIndex = 8;
            this.refresh_button.Text = "Refresh";
            this.refresh_button.UseVisualStyleBackColor = true;
            this.refresh_button.Click += new System.EventHandler(this.refresh_button_Click);
            // 
            // count_label
            // 
            this.count_label.AutoSize = true;
            this.count_label.Location = new System.Drawing.Point(353, 156);
            this.count_label.Name = "count_label";
            this.count_label.Size = new System.Drawing.Size(47, 13);
            this.count_label.TabIndex = 7;
            this.count_label.Text = "Count = ";
            // 
            // step_button
            // 
            this.step_button.Location = new System.Drawing.Point(341, 60);
            this.step_button.Name = "step_button";
            this.step_button.Size = new System.Drawing.Size(149, 83);
            this.step_button.TabIndex = 6;
            this.step_button.Text = "Step";
            this.step_button.UseVisualStyleBackColor = true;
            this.step_button.Click += new System.EventHandler(this.step_button_Click);
            // 
            // group1Dboundaries
            // 
            this.group1Dboundaries.Controls.Add(this.boundary1Ddescriptor_label);
            this.group1Dboundaries.Controls.Add(this.bottomV1Dval_label);
            this.group1Dboundaries.Controls.Add(this.topV1Dval);
            this.group1Dboundaries.Controls.Add(this.topV1Dval_label);
            this.group1Dboundaries.Controls.Add(this.bottomV1Dval);
            this.group1Dboundaries.Location = new System.Drawing.Point(35, 415);
            this.group1Dboundaries.Name = "group1Dboundaries";
            this.group1Dboundaries.Size = new System.Drawing.Size(263, 110);
            this.group1Dboundaries.TabIndex = 5;
            this.group1Dboundaries.TabStop = false;
            this.group1Dboundaries.Text = "Boundary Conditions";
            // 
            // boundary1Ddescriptor_label
            // 
            this.boundary1Ddescriptor_label.AutoSize = true;
            this.boundary1Ddescriptor_label.Location = new System.Drawing.Point(9, 75);
            this.boundary1Ddescriptor_label.Name = "boundary1Ddescriptor_label";
            this.boundary1Ddescriptor_label.Size = new System.Drawing.Size(226, 26);
            this.boundary1Ddescriptor_label.TabIndex = 16;
            this.boundary1Ddescriptor_label.Text = "Boundary conditions will be set to zero if empty\r\nThis should be the default valu" +
    "es";
            // 
            // bottomV1Dval_label
            // 
            this.bottomV1Dval_label.AutoSize = true;
            this.bottomV1Dval_label.Location = new System.Drawing.Point(118, 47);
            this.bottomV1Dval_label.Name = "bottomV1Dval_label";
            this.bottomV1Dval_label.Size = new System.Drawing.Size(139, 13);
            this.bottomV1Dval_label.TabIndex = 15;
            this.bottomV1Dval_label.Text = "Voltage at bottom of domain";
            // 
            // topV1Dval
            // 
            this.topV1Dval.Location = new System.Drawing.Point(6, 19);
            this.topV1Dval.Name = "topV1Dval";
            this.topV1Dval.Size = new System.Drawing.Size(94, 20);
            this.topV1Dval.TabIndex = 12;
            this.topV1Dval.Text = "0.0";
            // 
            // topV1Dval_label
            // 
            this.topV1Dval_label.AutoSize = true;
            this.topV1Dval_label.Location = new System.Drawing.Point(118, 21);
            this.topV1Dval_label.Name = "topV1Dval_label";
            this.topV1Dval_label.Size = new System.Drawing.Size(122, 13);
            this.topV1Dval_label.TabIndex = 13;
            this.topV1Dval_label.Text = "Voltage at top of domain";
            // 
            // bottomV1Dval
            // 
            this.bottomV1Dval.Location = new System.Drawing.Point(6, 45);
            this.bottomV1Dval.Name = "bottomV1Dval";
            this.bottomV1Dval.Size = new System.Drawing.Size(94, 20);
            this.bottomV1Dval.TabIndex = 14;
            this.bottomV1Dval.Text = "0.0";
            // 
            // groupConversion
            // 
            this.groupConversion.Controls.Add(this.ny1Dval_label);
            this.groupConversion.Controls.Add(this.nx1Dval_label);
            this.groupConversion.Controls.Add(this.ny1Dval);
            this.groupConversion.Controls.Add(this.nx1Dval);
            this.groupConversion.Location = new System.Drawing.Point(35, 325);
            this.groupConversion.Name = "groupConversion";
            this.groupConversion.Size = new System.Drawing.Size(263, 75);
            this.groupConversion.TabIndex = 4;
            this.groupConversion.TabStop = false;
            this.groupConversion.Text = "Dopent Conversion Parameters";
            // 
            // ny1Dval_label
            // 
            this.ny1Dval_label.AutoSize = true;
            this.ny1Dval_label.Location = new System.Drawing.Point(125, 47);
            this.ny1Dval_label.Name = "ny1Dval_label";
            this.ny1Dval_label.Size = new System.Drawing.Size(20, 13);
            this.ny1Dval_label.TabIndex = 11;
            this.ny1Dval_label.Text = "Ny";
            // 
            // nx1Dval_label
            // 
            this.nx1Dval_label.AutoSize = true;
            this.nx1Dval_label.Location = new System.Drawing.Point(125, 21);
            this.nx1Dval_label.Name = "nx1Dval_label";
            this.nx1Dval_label.Size = new System.Drawing.Size(20, 13);
            this.nx1Dval_label.TabIndex = 9;
            this.nx1Dval_label.Text = "Nx";
            // 
            // ny1Dval
            // 
            this.ny1Dval.Location = new System.Drawing.Point(13, 45);
            this.ny1Dval.Name = "ny1Dval";
            this.ny1Dval.Size = new System.Drawing.Size(94, 20);
            this.ny1Dval.TabIndex = 10;
            // 
            // nx1Dval
            // 
            this.nx1Dval.Location = new System.Drawing.Point(13, 19);
            this.nx1Dval.Name = "nx1Dval";
            this.nx1Dval.Size = new System.Drawing.Size(94, 20);
            this.nx1Dval.TabIndex = 8;
            // 
            // groupPhysical
            // 
            this.groupPhysical.Controls.Add(this.temperature1Dval_label);
            this.groupPhysical.Controls.Add(this.temperature1Dval);
            this.groupPhysical.Location = new System.Drawing.Point(36, 261);
            this.groupPhysical.Name = "groupPhysical";
            this.groupPhysical.Size = new System.Drawing.Size(262, 48);
            this.groupPhysical.TabIndex = 3;
            this.groupPhysical.TabStop = false;
            this.groupPhysical.Text = "Physical Parameters";
            // 
            // temperature1Dval_label
            // 
            this.temperature1Dval_label.AutoSize = true;
            this.temperature1Dval_label.Location = new System.Drawing.Point(124, 21);
            this.temperature1Dval_label.Name = "temperature1Dval_label";
            this.temperature1Dval_label.Size = new System.Drawing.Size(112, 13);
            this.temperature1Dval_label.TabIndex = 9;
            this.temperature1Dval_label.Text = "Base Temperature / K";
            // 
            // temperature1Dval
            // 
            this.temperature1Dval.Location = new System.Drawing.Point(12, 19);
            this.temperature1Dval.Name = "temperature1Dval";
            this.temperature1Dval.Size = new System.Drawing.Size(94, 20);
            this.temperature1Dval.TabIndex = 8;
            this.temperature1Dval.Text = "1.0";
            // 
            // groupDFT
            // 
            this.groupDFT.Controls.Add(this.dftzmin1Dval_label);
            this.groupDFT.Controls.Add(this.dftzmin1Dval);
            this.groupDFT.Controls.Add(this.dftnz1Dval_label);
            this.groupDFT.Controls.Add(this.dft1DCheck);
            this.groupDFT.Controls.Add(this.dftnz1Dval);
            this.groupDFT.Location = new System.Drawing.Point(36, 146);
            this.groupDFT.Name = "groupDFT";
            this.groupDFT.Size = new System.Drawing.Size(262, 102);
            this.groupDFT.TabIndex = 2;
            this.groupDFT.TabStop = false;
            this.groupDFT.Text = "DFT Calculation Parameters";
            // 
            // dftzmin1Dval_label
            // 
            this.dftzmin1Dval_label.AutoSize = true;
            this.dftzmin1Dval_label.Location = new System.Drawing.Point(124, 70);
            this.dftzmin1Dval_label.Name = "dftzmin1Dval_label";
            this.dftzmin1Dval_label.Size = new System.Drawing.Size(94, 13);
            this.dftzmin1Dval_label.TabIndex = 7;
            this.dftzmin1Dval_label.Text = "Zmin for DFT / nm";
            // 
            // dftzmin1Dval
            // 
            this.dftzmin1Dval.Location = new System.Drawing.Point(12, 68);
            this.dftzmin1Dval.Name = "dftzmin1Dval";
            this.dftzmin1Dval.Size = new System.Drawing.Size(94, 20);
            this.dftzmin1Dval.TabIndex = 6;
            this.dftzmin1Dval.Text = "-330.0";
            // 
            // dftnz1Dval_label
            // 
            this.dftnz1Dval_label.AutoSize = true;
            this.dftnz1Dval_label.Location = new System.Drawing.Point(124, 44);
            this.dftnz1Dval_label.Name = "dftnz1Dval_label";
            this.dftnz1Dval_label.Size = new System.Drawing.Size(59, 13);
            this.dftnz1Dval_label.TabIndex = 5;
            this.dftnz1Dval_label.Text = "Nz for DFT";
            // 
            // dft1DCheck
            // 
            this.dft1DCheck.AutoSize = true;
            this.dft1DCheck.Checked = true;
            this.dft1DCheck.CheckState = System.Windows.Forms.CheckState.Checked;
            this.dft1DCheck.Location = new System.Drawing.Point(12, 19);
            this.dft1DCheck.Name = "dft1DCheck";
            this.dft1DCheck.Size = new System.Drawing.Size(78, 17);
            this.dft1DCheck.TabIndex = 0;
            this.dft1DCheck.Text = "With DFT?";
            this.dft1DCheck.UseVisualStyleBackColor = true;
            this.dft1DCheck.CheckedChanged += new System.EventHandler(this.dftCheck_CheckedChanged);
            // 
            // dftnz1Dval
            // 
            this.dftnz1Dval.Location = new System.Drawing.Point(12, 42);
            this.dftnz1Dval.Name = "dftnz1Dval";
            this.dftnz1Dval.Size = new System.Drawing.Size(94, 20);
            this.dftnz1Dval.TabIndex = 4;
            this.dftnz1Dval.Text = "100";
            // 
            // groupPotential
            // 
            this.groupPotential.Controls.Add(this.nz1Dval_label);
            this.groupPotential.Controls.Add(this.nz1Dval);
            this.groupPotential.Controls.Add(this.dz1Dval_label);
            this.groupPotential.Controls.Add(this.dz1Dval);
            this.groupPotential.Location = new System.Drawing.Point(35, 60);
            this.groupPotential.Name = "groupPotential";
            this.groupPotential.Size = new System.Drawing.Size(263, 80);
            this.groupPotential.TabIndex = 1;
            this.groupPotential.TabStop = false;
            this.groupPotential.Text = "Potential Parameters";
            // 
            // nz1Dval_label
            // 
            this.nz1Dval_label.AutoSize = true;
            this.nz1Dval_label.Location = new System.Drawing.Point(125, 52);
            this.nz1Dval_label.Name = "nz1Dval_label";
            this.nz1Dval_label.Size = new System.Drawing.Size(20, 13);
            this.nz1Dval_label.TabIndex = 3;
            this.nz1Dval_label.Text = "Nz";
            // 
            // nz1Dval
            // 
            this.nz1Dval.Location = new System.Drawing.Point(13, 50);
            this.nz1Dval.Name = "nz1Dval";
            this.nz1Dval.Size = new System.Drawing.Size(94, 20);
            this.nz1Dval.TabIndex = 2;
            this.nz1Dval.Text = "1900";
            // 
            // dz1Dval_label
            // 
            this.dz1Dval_label.AutoSize = true;
            this.dz1Dval_label.Location = new System.Drawing.Point(125, 26);
            this.dz1Dval_label.Name = "dz1Dval_label";
            this.dz1Dval_label.Size = new System.Drawing.Size(43, 13);
            this.dz1Dval_label.TabIndex = 1;
            this.dz1Dval_label.Text = "dz / nm";
            // 
            // dz1Dval
            // 
            this.dz1Dval.Location = new System.Drawing.Point(13, 24);
            this.dz1Dval.Name = "dz1Dval";
            this.dz1Dval.Size = new System.Drawing.Size(94, 20);
            this.dz1Dval.TabIndex = 0;
            this.dz1Dval.Text = "1.0";
            // 
            // dopentDescriptor
            // 
            this.dopentDescriptor.AutoSize = true;
            this.dopentDescriptor.Location = new System.Drawing.Point(16, 15);
            this.dopentDescriptor.Name = "dopentDescriptor";
            this.dopentDescriptor.Size = new System.Drawing.Size(434, 26);
            this.dopentDescriptor.TabIndex = 0;
            this.dopentDescriptor.Text = "Parameters for running the initial 1D calculation in order to calculate the dopen" +
    "t distribution.\r\nNOTE: This is redundant for 1D simulations";
            // 
            // otherinputs
            // 
            this.otherinputs.Controls.Add(this.groupConvergence);
            this.otherinputs.Controls.Add(this.dimensionLabel);
            this.otherinputs.Controls.Add(this.dimensionality);
            this.otherinputs.Controls.Add(this.includebatchdata);
            this.otherinputs.Controls.Add(this.checkpointing);
            this.otherinputs.Location = new System.Drawing.Point(4, 22);
            this.otherinputs.Name = "otherinputs";
            this.otherinputs.Size = new System.Drawing.Size(721, 620);
            this.otherinputs.TabIndex = 4;
            this.otherinputs.Text = "Other";
            this.otherinputs.UseVisualStyleBackColor = true;
            // 
            // groupConvergence
            // 
            this.groupConvergence.Controls.Add(this.max_iterations1d_val);
            this.groupConvergence.Controls.Add(this.max_iterations1d_labal);
            this.groupConvergence.Controls.Add(this.dopent_convergenceparams_label);
            this.groupConvergence.Controls.Add(this.main_convergenceparams_label);
            this.groupConvergence.Controls.Add(this.maxiteration_val);
            this.groupConvergence.Controls.Add(this.max_iteration_label);
            this.groupConvergence.Controls.Add(this.toleranceval);
            this.groupConvergence.Controls.Add(this.tolerance_label);
            this.groupConvergence.Controls.Add(this.tolerance1Dval);
            this.groupConvergence.Controls.Add(this.tolerance1D_label);
            this.groupConvergence.Location = new System.Drawing.Point(215, 41);
            this.groupConvergence.Name = "groupConvergence";
            this.groupConvergence.Size = new System.Drawing.Size(480, 211);
            this.groupConvergence.TabIndex = 5;
            this.groupConvergence.TabStop = false;
            this.groupConvergence.Text = "Convergence Parameters";
            // 
            // max_iterations1d_val
            // 
            this.max_iterations1d_val.Location = new System.Drawing.Point(212, 77);
            this.max_iterations1d_val.Name = "max_iterations1d_val";
            this.max_iterations1d_val.Size = new System.Drawing.Size(96, 20);
            this.max_iterations1d_val.TabIndex = 13;
            this.max_iterations1d_val.Text = "1000";
            // 
            // max_iterations1d_labal
            // 
            this.max_iterations1d_labal.AutoSize = true;
            this.max_iterations1d_labal.Location = new System.Drawing.Point(314, 80);
            this.max_iterations1d_labal.Name = "max_iterations1d_labal";
            this.max_iterations1d_labal.Size = new System.Drawing.Size(73, 13);
            this.max_iterations1d_labal.TabIndex = 12;
            this.max_iterations1d_labal.Text = "Max Iterations";
            // 
            // dopent_convergenceparams_label
            // 
            this.dopent_convergenceparams_label.AutoSize = true;
            this.dopent_convergenceparams_label.Location = new System.Drawing.Point(209, 25);
            this.dopent_convergenceparams_label.Name = "dopent_convergenceparams_label";
            this.dopent_convergenceparams_label.Size = new System.Drawing.Size(93, 13);
            this.dopent_convergenceparams_label.TabIndex = 11;
            this.dopent_convergenceparams_label.Text = "Dopent Simulation";
            // 
            // main_convergenceparams_label
            // 
            this.main_convergenceparams_label.AutoSize = true;
            this.main_convergenceparams_label.Location = new System.Drawing.Point(14, 25);
            this.main_convergenceparams_label.Name = "main_convergenceparams_label";
            this.main_convergenceparams_label.Size = new System.Drawing.Size(81, 13);
            this.main_convergenceparams_label.TabIndex = 10;
            this.main_convergenceparams_label.Text = "Main Simulation";
            // 
            // maxiteration_val
            // 
            this.maxiteration_val.Location = new System.Drawing.Point(15, 76);
            this.maxiteration_val.Name = "maxiteration_val";
            this.maxiteration_val.Size = new System.Drawing.Size(96, 20);
            this.maxiteration_val.TabIndex = 9;
            this.maxiteration_val.Text = "1000";
            // 
            // max_iteration_label
            // 
            this.max_iteration_label.AutoSize = true;
            this.max_iteration_label.Location = new System.Drawing.Point(117, 79);
            this.max_iteration_label.Name = "max_iteration_label";
            this.max_iteration_label.Size = new System.Drawing.Size(73, 13);
            this.max_iteration_label.TabIndex = 8;
            this.max_iteration_label.Text = "Max Iterations";
            // 
            // toleranceval
            // 
            this.toleranceval.Location = new System.Drawing.Point(15, 50);
            this.toleranceval.Name = "toleranceval";
            this.toleranceval.Size = new System.Drawing.Size(96, 20);
            this.toleranceval.TabIndex = 7;
            this.toleranceval.Text = "1.0e-6";
            // 
            // tolerance_label
            // 
            this.tolerance_label.AutoSize = true;
            this.tolerance_label.Location = new System.Drawing.Point(117, 53);
            this.tolerance_label.Name = "tolerance_label";
            this.tolerance_label.Size = new System.Drawing.Size(55, 13);
            this.tolerance_label.TabIndex = 6;
            this.tolerance_label.Text = "Tolerance";
            // 
            // tolerance1Dval
            // 
            this.tolerance1Dval.Location = new System.Drawing.Point(212, 51);
            this.tolerance1Dval.Name = "tolerance1Dval";
            this.tolerance1Dval.Size = new System.Drawing.Size(96, 20);
            this.tolerance1Dval.TabIndex = 5;
            this.tolerance1Dval.Text = "1.0e-6";
            // 
            // tolerance1D_label
            // 
            this.tolerance1D_label.AutoSize = true;
            this.tolerance1D_label.Location = new System.Drawing.Point(314, 54);
            this.tolerance1D_label.Name = "tolerance1D_label";
            this.tolerance1D_label.Size = new System.Drawing.Size(55, 13);
            this.tolerance1D_label.TabIndex = 4;
            this.tolerance1D_label.Text = "Tolerance";
            // 
            // dimensionLabel
            // 
            this.dimensionLabel.AutoSize = true;
            this.dimensionLabel.Location = new System.Drawing.Point(17, 15);
            this.dimensionLabel.Name = "dimensionLabel";
            this.dimensionLabel.Size = new System.Drawing.Size(56, 13);
            this.dimensionLabel.TabIndex = 3;
            this.dimensionLabel.Text = "Dimension";
            // 
            // dimensionality
            // 
            this.dimensionality.FormattingEnabled = true;
            this.dimensionality.Items.AddRange(new object[] {
            "1D",
            "2D",
            "3D"});
            this.dimensionality.Location = new System.Drawing.Point(20, 31);
            this.dimensionality.Name = "dimensionality";
            this.dimensionality.Size = new System.Drawing.Size(120, 43);
            this.dimensionality.TabIndex = 2;
            this.dimensionality.SelectedIndexChanged += new System.EventHandler(this.dimensionality_SelectedIndexChanged);
            // 
            // includebatchdata
            // 
            this.includebatchdata.AutoSize = true;
            this.includebatchdata.Checked = true;
            this.includebatchdata.CheckState = System.Windows.Forms.CheckState.Checked;
            this.includebatchdata.Location = new System.Drawing.Point(20, 117);
            this.includebatchdata.Name = "includebatchdata";
            this.includebatchdata.Size = new System.Drawing.Size(118, 30);
            this.includebatchdata.TabIndex = 1;
            this.includebatchdata.Text = "Include Batch Data\r\nin Input File?";
            this.includebatchdata.UseVisualStyleBackColor = true;
            // 
            // checkpointing
            // 
            this.checkpointing.AutoSize = true;
            this.checkpointing.Location = new System.Drawing.Point(20, 94);
            this.checkpointing.Name = "checkpointing";
            this.checkpointing.Size = new System.Drawing.Size(128, 17);
            this.checkpointing.TabIndex = 0;
            this.checkpointing.Text = "Allow Checkpointing?";
            this.checkpointing.UseVisualStyleBackColor = true;
            // 
            // subsplit
            // 
            this.subsplit.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
            this.subsplit.Dock = System.Windows.Forms.DockStyle.Fill;
            this.subsplit.Location = new System.Drawing.Point(0, 0);
            this.subsplit.Name = "subsplit";
            this.subsplit.Orientation = System.Windows.Forms.Orientation.Horizontal;
            // 
            // subsplit.Panel1
            // 
            this.subsplit.Panel1.Controls.Add(this.conduction_band);
            this.subsplit.Panel1.RightToLeft = System.Windows.Forms.RightToLeft.No;
            // 
            // subsplit.Panel2
            // 
            this.subsplit.Panel2.Controls.Add(this.density);
            this.subsplit.Panel2.RightToLeft = System.Windows.Forms.RightToLeft.No;
            this.subsplit.RightToLeft = System.Windows.Forms.RightToLeft.No;
            this.subsplit.Size = new System.Drawing.Size(678, 649);
            this.subsplit.SplitterDistance = 332;
            this.subsplit.TabIndex = 0;
            // 
            // conduction_band
            // 
            chartArea1.Name = "ChartArea1";
            this.conduction_band.ChartAreas.Add(chartArea1);
            legend1.Name = "Legend1";
            this.conduction_band.Legends.Add(legend1);
            this.conduction_band.Location = new System.Drawing.Point(3, 3);
            this.conduction_band.Name = "conduction_band";
            series1.ChartArea = "ChartArea1";
            series1.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series1.Legend = "Legend1";
            series1.Name = "conduction_band_data";
            series2.ChartArea = "ChartArea1";
            series2.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series2.Legend = "Legend1";
            series2.Name = "valence_band_data";
            series3.ChartArea = "ChartArea1";
            series3.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series3.Legend = "Legend1";
            series3.Name = "x_data";
            this.conduction_band.Series.Add(series1);
            this.conduction_band.Series.Add(series2);
            this.conduction_band.Series.Add(series3);
            this.conduction_band.Size = new System.Drawing.Size(668, 322);
            this.conduction_band.TabIndex = 5;
            this.conduction_band.Text = "chart1";
            // 
            // density
            // 
            chartArea2.Name = "ChartArea1";
            this.density.ChartAreas.Add(chartArea2);
            legend2.Name = "Legend1";
            this.density.Legends.Add(legend2);
            this.density.Location = new System.Drawing.Point(3, 3);
            this.density.Name = "density";
            series4.ChartArea = "ChartArea1";
            series4.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series4.Legend = "Legend1";
            series4.Name = "car_dens_data";
            series5.ChartArea = "ChartArea1";
            series5.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series5.Legend = "Legend1";
            series5.Name = "dop_dens_data";
            series6.ChartArea = "ChartArea1";
            series6.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series6.Legend = "Legend1";
            series6.Name = "gphi_data";
            this.density.Series.Add(series4);
            this.density.Series.Add(series5);
            this.density.Series.Add(series6);
            this.density.Size = new System.Drawing.Size(667, 303);
            this.density.TabIndex = 0;
            this.density.Text = "chart1";
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // carrier_dopent_density_Text
            // 
            this.carrier_dopent_density_Text.Location = new System.Drawing.Point(519, 369);
            this.carrier_dopent_density_Text.Name = "carrier_dopent_density_Text";
            this.carrier_dopent_density_Text.ReadOnly = true;
            this.carrier_dopent_density_Text.Size = new System.Drawing.Size(103, 20);
            this.carrier_dopent_density_Text.TabIndex = 14;
            this.carrier_dopent_density_Text.Text = "0.0";
            this.carrier_dopent_density_Text.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // carrier_density_dopent_label
            // 
            this.carrier_density_dopent_label.AutoSize = true;
            this.carrier_density_dopent_label.Location = new System.Drawing.Point(304, 372);
            this.carrier_density_dopent_label.Name = "carrier_density_dopent_label";
            this.carrier_density_dopent_label.Size = new System.Drawing.Size(209, 13);
            this.carrier_density_dopent_label.TabIndex = 15;
            this.carrier_density_dopent_label.Text = "Carrier density at heterostructure interface: ";
            // 
            // carrier_density_dopent_units
            // 
            this.carrier_density_dopent_units.AutoSize = true;
            this.carrier_density_dopent_units.Location = new System.Drawing.Point(634, 373);
            this.carrier_density_dopent_units.Name = "carrier_density_dopent_units";
            this.carrier_density_dopent_units.Size = new System.Drawing.Size(36, 13);
            this.carrier_density_dopent_units.TabIndex = 16;
            this.carrier_density_dopent_units.Text = "cm^-2";
            // 
            // Solver_GUI
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1416, 673);
            this.Controls.Add(this.splitContainer);
            this.Controls.Add(this.menuStrip1);
            this.MainMenuStrip = this.menuStrip1;
            this.Name = "Solver_GUI";
            this.Text = "Solver - 2D";
            this.menuStrip1.ResumeLayout(false);
            this.menuStrip1.PerformLayout();
            this.splitContainer.Panel1.ResumeLayout(false);
            this.splitContainer.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer)).EndInit();
            this.splitContainer.ResumeLayout(false);
            this.inputTab.ResumeLayout(false);
            this.bandstructure.ResumeLayout(false);
            this.bandstructure.PerformLayout();
            this.group_newlayer.ResumeLayout(false);
            this.group_newlayer.PerformLayout();
            this.densityinputs.ResumeLayout(false);
            this.densityinputs.PerformLayout();
            this.dopentinputs.ResumeLayout(false);
            this.dopentinputs.PerformLayout();
            this.groupPlotting.ResumeLayout(false);
            this.groupPlotting.PerformLayout();
            this.group1Dboundaries.ResumeLayout(false);
            this.group1Dboundaries.PerformLayout();
            this.groupConversion.ResumeLayout(false);
            this.groupConversion.PerformLayout();
            this.groupPhysical.ResumeLayout(false);
            this.groupPhysical.PerformLayout();
            this.groupDFT.ResumeLayout(false);
            this.groupDFT.PerformLayout();
            this.groupPotential.ResumeLayout(false);
            this.groupPotential.PerformLayout();
            this.otherinputs.ResumeLayout(false);
            this.otherinputs.PerformLayout();
            this.groupConvergence.ResumeLayout(false);
            this.groupConvergence.PerformLayout();
            this.subsplit.Panel1.ResumeLayout(false);
            this.subsplit.Panel2.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.subsplit)).EndInit();
            this.subsplit.ResumeLayout(false);
            ((System.ComponentModel.ISupportInitialize)(this.conduction_band)).EndInit();
            ((System.ComponentModel.ISupportInitialize)(this.density)).EndInit();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.MenuStrip menuStrip1;
        private System.Windows.Forms.ToolStripMenuItem fileToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem mnuNew;
        private System.Windows.Forms.ToolStripMenuItem mnuNewBandStructure;
        private System.Windows.Forms.ToolStripMenuItem mnuNewInputFile;
        private System.Windows.Forms.ToolStripMenuItem mnuOpen;
        private System.Windows.Forms.ToolStripMenuItem mnuOpenBandStructure;
        private System.Windows.Forms.ToolStripMenuItem mnuOpenInputFile;
        private System.Windows.Forms.ToolStripMenuItem mnuExit;
        private System.Windows.Forms.ToolStripMenuItem viewToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem settingToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem helpToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem gettingStartedToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem aboutToolStripMenuItem;
        private System.Windows.Forms.ToolStripSeparator toolStripMenuItem1;
        private System.Windows.Forms.ToolStripMenuItem mnuSave;
        private System.Windows.Forms.ToolStripMenuItem mnuSaveAs;
        private System.Windows.Forms.ToolStripSeparator toolStripMenuItem2;
        private System.Windows.Forms.SplitContainer splitContainer;
        private System.Windows.Forms.TabControl inputTab;
        private System.Windows.Forms.TabPage bandstructure;
        private System.Windows.Forms.TabPage potentialinputs;
        private System.Windows.Forms.TabPage densityinputs;
        private System.Windows.Forms.TabPage dopentinputs;
        private System.Windows.Forms.TabPage otherinputs;
        private System.Windows.Forms.SplitContainer subsplit;
        private System.Windows.Forms.Label dimensionLabel;
        private System.Windows.Forms.ListBox dimensionality;
        private System.Windows.Forms.CheckBox includebatchdata;
        private System.Windows.Forms.CheckBox checkpointing;
        private System.Windows.Forms.Label dopentDescriptor;
        private System.Windows.Forms.GroupBox groupDFT;
        private System.Windows.Forms.GroupBox groupPotential;
        private System.Windows.Forms.Label nz1Dval_label;
        private System.Windows.Forms.TextBox nz1Dval;
        private System.Windows.Forms.Label dz1Dval_label;
        private System.Windows.Forms.TextBox dz1Dval;
        private System.Windows.Forms.CheckBox dft1DCheck;
        private System.Windows.Forms.GroupBox groupPhysical;
        private System.Windows.Forms.Label temperature1Dval_label;
        private System.Windows.Forms.TextBox temperature1Dval;
        private System.Windows.Forms.Label dftzmin1Dval_label;
        private System.Windows.Forms.TextBox dftzmin1Dval;
        private System.Windows.Forms.Label dftnz1Dval_label;
        private System.Windows.Forms.TextBox dftnz1Dval;
        private System.Windows.Forms.GroupBox groupConversion;
        private System.Windows.Forms.Label ny1Dval_label;
        private System.Windows.Forms.Label nx1Dval_label;
        private System.Windows.Forms.TextBox ny1Dval;
        private System.Windows.Forms.TextBox nx1Dval;
        private System.Windows.Forms.GroupBox group1Dboundaries;
        private System.Windows.Forms.Label boundary1Ddescriptor_label;
        private System.Windows.Forms.Label bottomV1Dval_label;
        private System.Windows.Forms.TextBox topV1Dval;
        private System.Windows.Forms.Label topV1Dval_label;
        private System.Windows.Forms.TextBox bottomV1Dval;
        private System.Windows.Forms.Button editLayer_Button;
        private System.Windows.Forms.Button deleteLayer_Button;
        private System.Windows.Forms.Button addlayer_Button;
        private System.Windows.Forms.Label bandstructureCombo_label;
        private System.Windows.Forms.ComboBox bandstructureCombo;
        private System.Windows.Forms.DataVisualization.Charting.Chart conduction_band;
        private System.Windows.Forms.DataVisualization.Charting.Chart density;
        private System.Windows.Forms.Label output_label;
        private System.Windows.Forms.Button step_button;
        private System.Windows.Forms.Button refresh_button;
        private System.Windows.Forms.Label count_label;
        private System.Windows.Forms.Label count_no_label;
        private System.Windows.Forms.Label temperature_val_label;
        private System.Windows.Forms.Label temperature_label;
        private System.Windows.Forms.Button run_button;
        private System.Windows.Forms.ToolStripMenuItem solverConfigToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem solverConfigToolStripMenuItem1;
        private System.Windows.Forms.GroupBox groupPlotting;
        private System.Windows.Forms.Label pot_ymax_label;
        private System.Windows.Forms.TextBox pot_ymax_val;
        private System.Windows.Forms.Label pot_ymin_label;
        private System.Windows.Forms.TextBox pot_ymin_val;
        private System.Windows.Forms.Label pot_xmax_label;
        private System.Windows.Forms.TextBox pot_xmax_val;
        private System.Windows.Forms.Label pot_plot_label;
        private System.Windows.Forms.Label pot_xmin_label;
        private System.Windows.Forms.TextBox pot_xmin_val;
        private System.Windows.Forms.Label dens_ymax_label;
        private System.Windows.Forms.TextBox dens_ymax_val;
        private System.Windows.Forms.Label dens_ymin_label;
        private System.Windows.Forms.TextBox dens_ymin_val;
        private System.Windows.Forms.Label dens_xmax_label;
        private System.Windows.Forms.TextBox dens_xmax_val;
        private System.Windows.Forms.Label dens_plot_label;
        private System.Windows.Forms.Label dens_xmin_label;
        private System.Windows.Forms.TextBox dens_xmin_val;
        private System.Windows.Forms.GroupBox groupConvergence;
        private System.Windows.Forms.TextBox toleranceval;
        private System.Windows.Forms.Label tolerance_label;
        private System.Windows.Forms.TextBox tolerance1Dval;
        private System.Windows.Forms.Label tolerance1D_label;
        private System.Windows.Forms.TextBox max_iterations1d_val;
        private System.Windows.Forms.Label max_iterations1d_labal;
        private System.Windows.Forms.Label dopent_convergenceparams_label;
        private System.Windows.Forms.Label main_convergenceparams_label;
        private System.Windows.Forms.TextBox maxiteration_val;
        private System.Windows.Forms.Label max_iteration_label;
        private System.Windows.Forms.ListView bandstructure_list;
        private System.Windows.Forms.ColumnHeader material_header;
        private System.Windows.Forms.ColumnHeader thickness_header;
        private System.Windows.Forms.ColumnHeader alloy_header;
        private System.Windows.Forms.ColumnHeader donor_header;
        private System.Windows.Forms.ColumnHeader acceptor_header;
        private System.Windows.Forms.GroupBox group_newlayer;
        private System.Windows.Forms.TextBox newlayer_naval;
        private System.Windows.Forms.Label newlayer_na_label;
        private System.Windows.Forms.TextBox newlayer_ndval;
        private System.Windows.Forms.Label newlayer_nd_label;
        private System.Windows.Forms.TextBox newlayer_xval;
        private System.Windows.Forms.Label newlayer_x_label;
        private System.Windows.Forms.TextBox newlayer_thicknessval;
        private System.Windows.Forms.Label newlayer_thickness_label;
        private System.Windows.Forms.Label material_label;
        private System.Windows.Forms.ComboBox material_combo;
        private System.Windows.Forms.Label layerselect_label;
        private System.Windows.Forms.TextBox bandstructurefilename;
        private System.Windows.Forms.Label bandstructurefile_label;
        private System.Windows.Forms.Button bandstructure_create_button;
        private System.Windows.Forms.OpenFileDialog openFileDialog1;
        private System.Windows.Forms.Label total_thickness_label;
        private System.Windows.Forms.TextBox total_thickness_val;
        private System.Windows.Forms.Label carrier_density_dopent_label;
        private System.Windows.Forms.TextBox carrier_dopent_density_Text;
        private System.Windows.Forms.Label carrier_density_dopent_units;
    }
}

