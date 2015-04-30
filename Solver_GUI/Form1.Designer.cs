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
            this.menuStrip1 = new System.Windows.Forms.MenuStrip();
            this.fileToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuNew = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuNewBandStructure = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuNewInputFile = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuOpen = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuOpenBandStructure = new System.Windows.Forms.ToolStripMenuItem();
            this.mnuOpenInputFile = new System.Windows.Forms.ToolStripMenuItem();
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
            this.potentialinputs = new System.Windows.Forms.TabPage();
            this.densityinputs = new System.Windows.Forms.TabPage();
            this.dopentinputs = new System.Windows.Forms.TabPage();
            this.otherinputs = new System.Windows.Forms.TabPage();
            this.dimensionLabel = new System.Windows.Forms.Label();
            this.dimensionality = new System.Windows.Forms.ListBox();
            this.includebatchdata = new System.Windows.Forms.CheckBox();
            this.checkpointing = new System.Windows.Forms.CheckBox();
            this.subsplit = new System.Windows.Forms.SplitContainer();
            this.dopentDescriptor = new System.Windows.Forms.Label();
            this.groupPotential = new System.Windows.Forms.GroupBox();
            this.groupDFT = new System.Windows.Forms.GroupBox();
            this.dz1Dval = new System.Windows.Forms.TextBox();
            this.dz1Dval_label = new System.Windows.Forms.Label();
            this.nz1Dval_label = new System.Windows.Forms.Label();
            this.nz1Dval = new System.Windows.Forms.TextBox();
            this.dftCheck = new System.Windows.Forms.CheckBox();
            this.dftnz1Dval_label = new System.Windows.Forms.Label();
            this.dftnz1Dval = new System.Windows.Forms.TextBox();
            this.dftzmin1Dval_label = new System.Windows.Forms.Label();
            this.dftzmin1Dval = new System.Windows.Forms.TextBox();
            this.groupPhysical = new System.Windows.Forms.GroupBox();
            this.temperature1Dval_label = new System.Windows.Forms.Label();
            this.temperature1Dval = new System.Windows.Forms.TextBox();
            this.groupConversion = new System.Windows.Forms.GroupBox();
            this.ny1Dval_label = new System.Windows.Forms.Label();
            this.ny1Dval = new System.Windows.Forms.TextBox();
            this.nx1Dval_label = new System.Windows.Forms.Label();
            this.nx1Dval = new System.Windows.Forms.TextBox();
            this.group1Dboundaries = new System.Windows.Forms.GroupBox();
            this.bottomV1Dval_label = new System.Windows.Forms.Label();
            this.topV1Dval_label = new System.Windows.Forms.Label();
            this.bottomV1Dval = new System.Windows.Forms.TextBox();
            this.topV1Dval = new System.Windows.Forms.TextBox();
            this.boundary1Ddescriptor_label = new System.Windows.Forms.Label();
            this.bandstructureCombo = new System.Windows.Forms.ComboBox();
            this.bandstructureCombo_label = new System.Windows.Forms.Label();
            this.addlayer_Button = new System.Windows.Forms.Button();
            this.deleteLayer_Button = new System.Windows.Forms.Button();
            this.editLayer_Button = new System.Windows.Forms.Button();
            this.menuStrip1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.splitContainer)).BeginInit();
            this.splitContainer.Panel1.SuspendLayout();
            this.splitContainer.Panel2.SuspendLayout();
            this.splitContainer.SuspendLayout();
            this.inputTab.SuspendLayout();
            this.bandstructure.SuspendLayout();
            this.dopentinputs.SuspendLayout();
            this.otherinputs.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.subsplit)).BeginInit();
            this.subsplit.SuspendLayout();
            this.groupPotential.SuspendLayout();
            this.groupDFT.SuspendLayout();
            this.groupPhysical.SuspendLayout();
            this.groupConversion.SuspendLayout();
            this.group1Dboundaries.SuspendLayout();
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
            this.mnuNewInputFile});
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
            // mnuOpen
            // 
            this.mnuOpen.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.mnuOpenBandStructure,
            this.mnuOpenInputFile});
            this.mnuOpen.Name = "mnuOpen";
            this.mnuOpen.Size = new System.Drawing.Size(138, 22);
            this.mnuOpen.Text = "Open";
            // 
            // mnuOpenBandStructure
            // 
            this.mnuOpenBandStructure.Name = "mnuOpenBandStructure";
            this.mnuOpenBandStructure.Size = new System.Drawing.Size(152, 22);
            this.mnuOpenBandStructure.Text = "Band Structure";
            // 
            // mnuOpenInputFile
            // 
            this.mnuOpenInputFile.Name = "mnuOpenInputFile";
            this.mnuOpenInputFile.Size = new System.Drawing.Size(152, 22);
            this.mnuOpenInputFile.Text = "Input File";
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
            this.densityinputs.Location = new System.Drawing.Point(4, 22);
            this.densityinputs.Name = "densityinputs";
            this.densityinputs.Size = new System.Drawing.Size(721, 620);
            this.densityinputs.TabIndex = 2;
            this.densityinputs.Text = "Density";
            this.densityinputs.UseVisualStyleBackColor = true;
            // 
            // dopentinputs
            // 
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
            // otherinputs
            // 
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
            this.subsplit.Panel1.RightToLeft = System.Windows.Forms.RightToLeft.No;
            // 
            // subsplit.Panel2
            // 
            this.subsplit.Panel2.RightToLeft = System.Windows.Forms.RightToLeft.No;
            this.subsplit.RightToLeft = System.Windows.Forms.RightToLeft.No;
            this.subsplit.Size = new System.Drawing.Size(678, 649);
            this.subsplit.SplitterDistance = 332;
            this.subsplit.TabIndex = 0;
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
            // groupPotential
            // 
            this.groupPotential.Controls.Add(this.nz1Dval_label);
            this.groupPotential.Controls.Add(this.nz1Dval);
            this.groupPotential.Controls.Add(this.dz1Dval_label);
            this.groupPotential.Controls.Add(this.dz1Dval);
            this.groupPotential.Location = new System.Drawing.Point(40, 74);
            this.groupPotential.Name = "groupPotential";
            this.groupPotential.Size = new System.Drawing.Size(263, 80);
            this.groupPotential.TabIndex = 1;
            this.groupPotential.TabStop = false;
            this.groupPotential.Text = "Potential Parameters";
            // 
            // groupDFT
            // 
            this.groupDFT.Controls.Add(this.dftzmin1Dval_label);
            this.groupDFT.Controls.Add(this.dftzmin1Dval);
            this.groupDFT.Controls.Add(this.dftnz1Dval_label);
            this.groupDFT.Controls.Add(this.dftCheck);
            this.groupDFT.Controls.Add(this.dftnz1Dval);
            this.groupDFT.Location = new System.Drawing.Point(41, 160);
            this.groupDFT.Name = "groupDFT";
            this.groupDFT.Size = new System.Drawing.Size(262, 102);
            this.groupDFT.TabIndex = 2;
            this.groupDFT.TabStop = false;
            this.groupDFT.Text = "DFT Calculation Parameters";
            // 
            // dz1Dval
            // 
            this.dz1Dval.Location = new System.Drawing.Point(13, 24);
            this.dz1Dval.Name = "dz1Dval";
            this.dz1Dval.Size = new System.Drawing.Size(94, 20);
            this.dz1Dval.TabIndex = 0;
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
            // 
            // dftCheck
            // 
            this.dftCheck.AutoSize = true;
            this.dftCheck.Checked = true;
            this.dftCheck.CheckState = System.Windows.Forms.CheckState.Checked;
            this.dftCheck.Location = new System.Drawing.Point(12, 19);
            this.dftCheck.Name = "dftCheck";
            this.dftCheck.Size = new System.Drawing.Size(78, 17);
            this.dftCheck.TabIndex = 0;
            this.dftCheck.Text = "With DFT?";
            this.dftCheck.UseVisualStyleBackColor = true;
            this.dftCheck.CheckedChanged += new System.EventHandler(this.dftCheck_CheckedChanged);
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
            // dftnz1Dval
            // 
            this.dftnz1Dval.Location = new System.Drawing.Point(12, 42);
            this.dftnz1Dval.Name = "dftnz1Dval";
            this.dftnz1Dval.Size = new System.Drawing.Size(94, 20);
            this.dftnz1Dval.TabIndex = 4;
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
            // 
            // groupPhysical
            // 
            this.groupPhysical.Controls.Add(this.temperature1Dval_label);
            this.groupPhysical.Controls.Add(this.temperature1Dval);
            this.groupPhysical.Location = new System.Drawing.Point(41, 275);
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
            // groupConversion
            // 
            this.groupConversion.Controls.Add(this.ny1Dval_label);
            this.groupConversion.Controls.Add(this.nx1Dval_label);
            this.groupConversion.Controls.Add(this.ny1Dval);
            this.groupConversion.Controls.Add(this.nx1Dval);
            this.groupConversion.Location = new System.Drawing.Point(40, 339);
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
            // ny1Dval
            // 
            this.ny1Dval.Location = new System.Drawing.Point(13, 45);
            this.ny1Dval.Name = "ny1Dval";
            this.ny1Dval.Size = new System.Drawing.Size(94, 20);
            this.ny1Dval.TabIndex = 10;
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
            // nx1Dval
            // 
            this.nx1Dval.Location = new System.Drawing.Point(13, 19);
            this.nx1Dval.Name = "nx1Dval";
            this.nx1Dval.Size = new System.Drawing.Size(94, 20);
            this.nx1Dval.TabIndex = 8;
            // 
            // group1Dboundaries
            // 
            this.group1Dboundaries.Controls.Add(this.boundary1Ddescriptor_label);
            this.group1Dboundaries.Controls.Add(this.bottomV1Dval_label);
            this.group1Dboundaries.Controls.Add(this.topV1Dval);
            this.group1Dboundaries.Controls.Add(this.topV1Dval_label);
            this.group1Dboundaries.Controls.Add(this.bottomV1Dval);
            this.group1Dboundaries.Location = new System.Drawing.Point(40, 429);
            this.group1Dboundaries.Name = "group1Dboundaries";
            this.group1Dboundaries.Size = new System.Drawing.Size(263, 110);
            this.group1Dboundaries.TabIndex = 5;
            this.group1Dboundaries.TabStop = false;
            this.group1Dboundaries.Text = "Boundary Conditions";
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
            // 
            // topV1Dval
            // 
            this.topV1Dval.Location = new System.Drawing.Point(6, 19);
            this.topV1Dval.Name = "topV1Dval";
            this.topV1Dval.Size = new System.Drawing.Size(94, 20);
            this.topV1Dval.TabIndex = 12;
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
            // bandstructureCombo_label
            // 
            this.bandstructureCombo_label.AutoSize = true;
            this.bandstructureCombo_label.Location = new System.Drawing.Point(18, 22);
            this.bandstructureCombo_label.Name = "bandstructureCombo_label";
            this.bandstructureCombo_label.Size = new System.Drawing.Size(115, 13);
            this.bandstructureCombo_label.TabIndex = 1;
            this.bandstructureCombo_label.Text = "Current Band Structure";
            // 
            // addlayer_Button
            // 
            this.addlayer_Button.Location = new System.Drawing.Point(148, 37);
            this.addlayer_Button.Name = "addlayer_Button";
            this.addlayer_Button.Size = new System.Drawing.Size(91, 22);
            this.addlayer_Button.TabIndex = 2;
            this.addlayer_Button.Text = "Add Layer";
            this.addlayer_Button.UseVisualStyleBackColor = true;
            // 
            // deleteLayer_Button
            // 
            this.deleteLayer_Button.Location = new System.Drawing.Point(245, 37);
            this.deleteLayer_Button.Name = "deleteLayer_Button";
            this.deleteLayer_Button.Size = new System.Drawing.Size(91, 22);
            this.deleteLayer_Button.TabIndex = 3;
            this.deleteLayer_Button.Text = "Delete Layer";
            this.deleteLayer_Button.UseVisualStyleBackColor = true;
            // 
            // editLayer_Button
            // 
            this.editLayer_Button.Location = new System.Drawing.Point(342, 37);
            this.editLayer_Button.Name = "editLayer_Button";
            this.editLayer_Button.Size = new System.Drawing.Size(91, 22);
            this.editLayer_Button.TabIndex = 4;
            this.editLayer_Button.Text = "Edit Layer";
            this.editLayer_Button.UseVisualStyleBackColor = true;
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
            this.dopentinputs.ResumeLayout(false);
            this.dopentinputs.PerformLayout();
            this.otherinputs.ResumeLayout(false);
            this.otherinputs.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.subsplit)).EndInit();
            this.subsplit.ResumeLayout(false);
            this.groupPotential.ResumeLayout(false);
            this.groupPotential.PerformLayout();
            this.groupDFT.ResumeLayout(false);
            this.groupDFT.PerformLayout();
            this.groupPhysical.ResumeLayout(false);
            this.groupPhysical.PerformLayout();
            this.groupConversion.ResumeLayout(false);
            this.groupConversion.PerformLayout();
            this.group1Dboundaries.ResumeLayout(false);
            this.group1Dboundaries.PerformLayout();
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
        private System.Windows.Forms.CheckBox dftCheck;
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
    }
}

