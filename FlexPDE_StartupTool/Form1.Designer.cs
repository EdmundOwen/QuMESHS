namespace FlexPDE_StartupTool
{
    partial class Form1
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
            this.Start_button = new System.Windows.Forms.Button();
            this.Stop_button = new System.Windows.Forms.Button();
            this.timetostart_progressbar = new System.Windows.Forms.ProgressBar();
            this.time_text = new System.Windows.Forms.TextBox();
            this.Time_label = new System.Windows.Forms.Label();
            this.s_label = new System.Windows.Forms.Label();
            this.error_label = new System.Windows.Forms.Label();
            this.flexpde_filename = new System.Windows.Forms.TextBox();
            this.flexPDE_location = new System.Windows.Forms.Label();
            this.executable_select_button = new System.Windows.Forms.Button();
            this.openFileDialog1 = new System.Windows.Forms.OpenFileDialog();
            this.script_file_label = new System.Windows.Forms.Label();
            this.launch_name_1 = new System.Windows.Forms.TextBox();
            this.launch_name_2 = new System.Windows.Forms.TextBox();
            this.launch_name_3 = new System.Windows.Forms.TextBox();
            this.launch_name_4 = new System.Windows.Forms.TextBox();
            this.launch_name_5 = new System.Windows.Forms.TextBox();
            this.launch_check_1 = new System.Windows.Forms.CheckBox();
            this.launch_check_2 = new System.Windows.Forms.CheckBox();
            this.launch_check_3 = new System.Windows.Forms.CheckBox();
            this.launch_check_4 = new System.Windows.Forms.CheckBox();
            this.launch_check_5 = new System.Windows.Forms.CheckBox();
            this.launch_select_1 = new System.Windows.Forms.Button();
            this.launch_select_2 = new System.Windows.Forms.Button();
            this.launch_select_3 = new System.Windows.Forms.Button();
            this.launch_select_4 = new System.Windows.Forms.Button();
            this.launch_select_5 = new System.Windows.Forms.Button();
            this.run_check = new System.Windows.Forms.CheckBox();
            this.exit_check = new System.Windows.Forms.CheckBox();
            this.minimized_check = new System.Windows.Forms.CheckBox();
            this.silent_check = new System.Windows.Forms.CheckBox();
            this.launch_select_6 = new System.Windows.Forms.Button();
            this.launch_check_6 = new System.Windows.Forms.CheckBox();
            this.launch_name_6 = new System.Windows.Forms.TextBox();
            this.backgroundworker = new System.ComponentModel.BackgroundWorker();
            this.logout_button = new System.Windows.Forms.Button();
            this.session_no_text = new System.Windows.Forms.TextBox();
            this.session_id_label = new System.Windows.Forms.Label();
            this.SuspendLayout();
            // 
            // Start_button
            // 
            this.Start_button.Location = new System.Drawing.Point(12, 371);
            this.Start_button.Name = "Start_button";
            this.Start_button.Size = new System.Drawing.Size(75, 23);
            this.Start_button.TabIndex = 0;
            this.Start_button.Text = "Start";
            this.Start_button.UseVisualStyleBackColor = true;
            this.Start_button.Click += new System.EventHandler(this.start_Click);
            // 
            // Stop_button
            // 
            this.Stop_button.Location = new System.Drawing.Point(93, 371);
            this.Stop_button.Name = "Stop_button";
            this.Stop_button.Size = new System.Drawing.Size(75, 23);
            this.Stop_button.TabIndex = 1;
            this.Stop_button.Text = "Stop";
            this.Stop_button.UseVisualStyleBackColor = true;
            this.Stop_button.Click += new System.EventHandler(this.stop_Click);
            // 
            // timetostart_progressbar
            // 
            this.timetostart_progressbar.Location = new System.Drawing.Point(12, 341);
            this.timetostart_progressbar.Name = "timetostart_progressbar";
            this.timetostart_progressbar.Size = new System.Drawing.Size(382, 24);
            this.timetostart_progressbar.TabIndex = 2;
            // 
            // time_text
            // 
            this.time_text.Location = new System.Drawing.Point(242, 315);
            this.time_text.Name = "time_text";
            this.time_text.Size = new System.Drawing.Size(134, 20);
            this.time_text.TabIndex = 3;
            this.time_text.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // Time_label
            // 
            this.Time_label.AutoSize = true;
            this.Time_label.Location = new System.Drawing.Point(169, 318);
            this.Time_label.Name = "Time_label";
            this.Time_label.Size = new System.Drawing.Size(67, 13);
            this.Time_label.TabIndex = 4;
            this.Time_label.Text = "Time to Start";
            // 
            // s_label
            // 
            this.s_label.AutoSize = true;
            this.s_label.Location = new System.Drawing.Point(382, 318);
            this.s_label.Name = "s_label";
            this.s_label.Size = new System.Drawing.Size(12, 13);
            this.s_label.TabIndex = 5;
            this.s_label.Text = "s";
            // 
            // error_label
            // 
            this.error_label.AutoSize = true;
            this.error_label.Location = new System.Drawing.Point(93, 376);
            this.error_label.Name = "error_label";
            this.error_label.Size = new System.Drawing.Size(0, 13);
            this.error_label.TabIndex = 6;
            // 
            // flexpde_filename
            // 
            this.flexpde_filename.Location = new System.Drawing.Point(12, 25);
            this.flexpde_filename.Name = "flexpde_filename";
            this.flexpde_filename.Size = new System.Drawing.Size(268, 20);
            this.flexpde_filename.TabIndex = 7;
            this.flexpde_filename.Text = "C:\\FlexPDE6\\FlexPDE6.exe";
            // 
            // flexPDE_location
            // 
            this.flexPDE_location.AutoSize = true;
            this.flexPDE_location.Location = new System.Drawing.Point(12, 9);
            this.flexPDE_location.Name = "flexPDE_location";
            this.flexPDE_location.Size = new System.Drawing.Size(143, 13);
            this.flexPDE_location.TabIndex = 8;
            this.flexPDE_location.Text = "FlexPDE executable location";
            // 
            // executable_select_button
            // 
            this.executable_select_button.Location = new System.Drawing.Point(286, 22);
            this.executable_select_button.Name = "executable_select_button";
            this.executable_select_button.Size = new System.Drawing.Size(65, 23);
            this.executable_select_button.TabIndex = 9;
            this.executable_select_button.Text = "Select";
            this.executable_select_button.UseVisualStyleBackColor = true;
            this.executable_select_button.Click += new System.EventHandler(this.select_Click);
            // 
            // openFileDialog1
            // 
            this.openFileDialog1.FileName = "openFileDialog1";
            // 
            // script_file_label
            // 
            this.script_file_label.AutoSize = true;
            this.script_file_label.Location = new System.Drawing.Point(12, 60);
            this.script_file_label.Name = "script_file_label";
            this.script_file_label.Size = new System.Drawing.Size(58, 13);
            this.script_file_label.TabIndex = 10;
            this.script_file_label.Text = "Script Files";
            // 
            // launch_name_1
            // 
            this.launch_name_1.Location = new System.Drawing.Point(8, 76);
            this.launch_name_1.Name = "launch_name_1";
            this.launch_name_1.Size = new System.Drawing.Size(241, 20);
            this.launch_name_1.TabIndex = 11;
            // 
            // launch_name_2
            // 
            this.launch_name_2.Location = new System.Drawing.Point(8, 102);
            this.launch_name_2.Name = "launch_name_2";
            this.launch_name_2.Size = new System.Drawing.Size(242, 20);
            this.launch_name_2.TabIndex = 12;
            // 
            // launch_name_3
            // 
            this.launch_name_3.Location = new System.Drawing.Point(8, 128);
            this.launch_name_3.Name = "launch_name_3";
            this.launch_name_3.Size = new System.Drawing.Size(241, 20);
            this.launch_name_3.TabIndex = 13;
            // 
            // launch_name_4
            // 
            this.launch_name_4.Location = new System.Drawing.Point(8, 154);
            this.launch_name_4.Name = "launch_name_4";
            this.launch_name_4.Size = new System.Drawing.Size(242, 20);
            this.launch_name_4.TabIndex = 14;
            // 
            // launch_name_5
            // 
            this.launch_name_5.Location = new System.Drawing.Point(8, 180);
            this.launch_name_5.Name = "launch_name_5";
            this.launch_name_5.Size = new System.Drawing.Size(241, 20);
            this.launch_name_5.TabIndex = 15;
            // 
            // launch_check_1
            // 
            this.launch_check_1.AutoSize = true;
            this.launch_check_1.Checked = true;
            this.launch_check_1.CheckState = System.Windows.Forms.CheckState.Checked;
            this.launch_check_1.Location = new System.Drawing.Point(326, 78);
            this.launch_check_1.Name = "launch_check_1";
            this.launch_check_1.Size = new System.Drawing.Size(68, 17);
            this.launch_check_1.TabIndex = 16;
            this.launch_check_1.Text = "Launch?";
            this.launch_check_1.UseVisualStyleBackColor = true;
            // 
            // launch_check_2
            // 
            this.launch_check_2.AutoSize = true;
            this.launch_check_2.Location = new System.Drawing.Point(326, 104);
            this.launch_check_2.Name = "launch_check_2";
            this.launch_check_2.Size = new System.Drawing.Size(68, 17);
            this.launch_check_2.TabIndex = 17;
            this.launch_check_2.Text = "Launch?";
            this.launch_check_2.UseVisualStyleBackColor = true;
            // 
            // launch_check_3
            // 
            this.launch_check_3.AutoSize = true;
            this.launch_check_3.Location = new System.Drawing.Point(326, 130);
            this.launch_check_3.Name = "launch_check_3";
            this.launch_check_3.Size = new System.Drawing.Size(68, 17);
            this.launch_check_3.TabIndex = 18;
            this.launch_check_3.Text = "Launch?";
            this.launch_check_3.UseVisualStyleBackColor = true;
            // 
            // launch_check_4
            // 
            this.launch_check_4.AutoSize = true;
            this.launch_check_4.Location = new System.Drawing.Point(326, 156);
            this.launch_check_4.Name = "launch_check_4";
            this.launch_check_4.Size = new System.Drawing.Size(68, 17);
            this.launch_check_4.TabIndex = 19;
            this.launch_check_4.Text = "Launch?";
            this.launch_check_4.UseVisualStyleBackColor = true;
            // 
            // launch_check_5
            // 
            this.launch_check_5.AutoSize = true;
            this.launch_check_5.Location = new System.Drawing.Point(326, 182);
            this.launch_check_5.Name = "launch_check_5";
            this.launch_check_5.Size = new System.Drawing.Size(68, 17);
            this.launch_check_5.TabIndex = 20;
            this.launch_check_5.Text = "Launch?";
            this.launch_check_5.UseVisualStyleBackColor = true;
            // 
            // launch_select_1
            // 
            this.launch_select_1.Location = new System.Drawing.Point(255, 74);
            this.launch_select_1.Name = "launch_select_1";
            this.launch_select_1.Size = new System.Drawing.Size(65, 23);
            this.launch_select_1.TabIndex = 21;
            this.launch_select_1.Text = "Select";
            this.launch_select_1.UseVisualStyleBackColor = true;
            this.launch_select_1.Click += new System.EventHandler(this.launch_select_1_Click);
            // 
            // launch_select_2
            // 
            this.launch_select_2.Location = new System.Drawing.Point(255, 100);
            this.launch_select_2.Name = "launch_select_2";
            this.launch_select_2.Size = new System.Drawing.Size(65, 23);
            this.launch_select_2.TabIndex = 22;
            this.launch_select_2.Text = "Select";
            this.launch_select_2.UseVisualStyleBackColor = true;
            this.launch_select_2.Click += new System.EventHandler(this.launch_select_2_Click);
            // 
            // launch_select_3
            // 
            this.launch_select_3.Location = new System.Drawing.Point(255, 126);
            this.launch_select_3.Name = "launch_select_3";
            this.launch_select_3.Size = new System.Drawing.Size(65, 23);
            this.launch_select_3.TabIndex = 23;
            this.launch_select_3.Text = "Select";
            this.launch_select_3.UseVisualStyleBackColor = true;
            this.launch_select_3.Click += new System.EventHandler(this.launch_select_3_Click);
            // 
            // launch_select_4
            // 
            this.launch_select_4.Location = new System.Drawing.Point(255, 152);
            this.launch_select_4.Name = "launch_select_4";
            this.launch_select_4.Size = new System.Drawing.Size(65, 23);
            this.launch_select_4.TabIndex = 24;
            this.launch_select_4.Text = "Select";
            this.launch_select_4.UseVisualStyleBackColor = true;
            this.launch_select_4.Click += new System.EventHandler(this.launch_select_4_Click);
            // 
            // launch_select_5
            // 
            this.launch_select_5.Location = new System.Drawing.Point(255, 178);
            this.launch_select_5.Name = "launch_select_5";
            this.launch_select_5.Size = new System.Drawing.Size(65, 23);
            this.launch_select_5.TabIndex = 25;
            this.launch_select_5.Text = "Select";
            this.launch_select_5.UseVisualStyleBackColor = true;
            this.launch_select_5.Click += new System.EventHandler(this.launch_select_5_Click);
            // 
            // run_check
            // 
            this.run_check.AutoSize = true;
            this.run_check.Location = new System.Drawing.Point(8, 249);
            this.run_check.Name = "run_check";
            this.run_check.Size = new System.Drawing.Size(60, 17);
            this.run_check.TabIndex = 26;
            this.run_check.Text = "Run -R";
            this.run_check.UseVisualStyleBackColor = true;
            // 
            // exit_check
            // 
            this.exit_check.AutoSize = true;
            this.exit_check.Location = new System.Drawing.Point(8, 272);
            this.exit_check.Name = "exit_check";
            this.exit_check.Size = new System.Drawing.Size(56, 17);
            this.exit_check.TabIndex = 27;
            this.exit_check.Text = "Exit -X";
            this.exit_check.UseVisualStyleBackColor = true;
            // 
            // minimized_check
            // 
            this.minimized_check.AutoSize = true;
            this.minimized_check.Location = new System.Drawing.Point(8, 295);
            this.minimized_check.Name = "minimized_check";
            this.minimized_check.Size = new System.Drawing.Size(87, 17);
            this.minimized_check.TabIndex = 28;
            this.minimized_check.Text = "Minimized -M";
            this.minimized_check.UseVisualStyleBackColor = true;
            // 
            // silent_check
            // 
            this.silent_check.AutoSize = true;
            this.silent_check.Checked = true;
            this.silent_check.CheckState = System.Windows.Forms.CheckState.Checked;
            this.silent_check.Location = new System.Drawing.Point(8, 318);
            this.silent_check.Name = "silent_check";
            this.silent_check.Size = new System.Drawing.Size(65, 17);
            this.silent_check.TabIndex = 29;
            this.silent_check.Text = "Silent -S";
            this.silent_check.TextImageRelation = System.Windows.Forms.TextImageRelation.ImageBeforeText;
            this.silent_check.UseVisualStyleBackColor = true;
            // 
            // launch_select_6
            // 
            this.launch_select_6.Location = new System.Drawing.Point(255, 204);
            this.launch_select_6.Name = "launch_select_6";
            this.launch_select_6.Size = new System.Drawing.Size(65, 23);
            this.launch_select_6.TabIndex = 32;
            this.launch_select_6.Text = "Select";
            this.launch_select_6.UseVisualStyleBackColor = true;
            this.launch_select_6.Click += new System.EventHandler(this.launch_select_6_Click);
            // 
            // launch_check_6
            // 
            this.launch_check_6.AutoSize = true;
            this.launch_check_6.Location = new System.Drawing.Point(326, 208);
            this.launch_check_6.Name = "launch_check_6";
            this.launch_check_6.Size = new System.Drawing.Size(68, 17);
            this.launch_check_6.TabIndex = 31;
            this.launch_check_6.Text = "Launch?";
            this.launch_check_6.UseVisualStyleBackColor = true;
            // 
            // launch_name_6
            // 
            this.launch_name_6.Location = new System.Drawing.Point(8, 206);
            this.launch_name_6.Name = "launch_name_6";
            this.launch_name_6.Size = new System.Drawing.Size(241, 20);
            this.launch_name_6.TabIndex = 30;
            // 
            // backgroundworker
            // 
            this.backgroundworker.WorkerReportsProgress = true;
            this.backgroundworker.WorkerSupportsCancellation = true;
            // 
            // logout_button
            // 
            this.logout_button.Location = new System.Drawing.Point(315, 371);
            this.logout_button.Name = "logout_button";
            this.logout_button.Size = new System.Drawing.Size(75, 23);
            this.logout_button.TabIndex = 33;
            this.logout_button.Text = "Log Out";
            this.logout_button.UseVisualStyleBackColor = true;
            this.logout_button.Click += new System.EventHandler(this.logout_button_Click);
            // 
            // session_no_text
            // 
            this.session_no_text.Location = new System.Drawing.Point(286, 373);
            this.session_no_text.Name = "session_no_text";
            this.session_no_text.Size = new System.Drawing.Size(23, 20);
            this.session_no_text.TabIndex = 34;
            this.session_no_text.Text = "1";
            this.session_no_text.TextAlign = System.Windows.Forms.HorizontalAlignment.Right;
            // 
            // session_id_label
            // 
            this.session_id_label.AutoSize = true;
            this.session_id_label.Location = new System.Drawing.Point(219, 376);
            this.session_id_label.Name = "session_id_label";
            this.session_id_label.Size = new System.Drawing.Size(61, 13);
            this.session_id_label.TabIndex = 35;
            this.session_id_label.Text = "Session ID:";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(402, 402);
            this.Controls.Add(this.session_id_label);
            this.Controls.Add(this.session_no_text);
            this.Controls.Add(this.logout_button);
            this.Controls.Add(this.launch_select_6);
            this.Controls.Add(this.launch_check_6);
            this.Controls.Add(this.launch_name_6);
            this.Controls.Add(this.silent_check);
            this.Controls.Add(this.minimized_check);
            this.Controls.Add(this.exit_check);
            this.Controls.Add(this.run_check);
            this.Controls.Add(this.launch_select_5);
            this.Controls.Add(this.launch_select_4);
            this.Controls.Add(this.launch_select_3);
            this.Controls.Add(this.launch_select_2);
            this.Controls.Add(this.launch_select_1);
            this.Controls.Add(this.launch_check_5);
            this.Controls.Add(this.launch_check_4);
            this.Controls.Add(this.launch_check_3);
            this.Controls.Add(this.launch_check_2);
            this.Controls.Add(this.launch_check_1);
            this.Controls.Add(this.launch_name_5);
            this.Controls.Add(this.launch_name_4);
            this.Controls.Add(this.launch_name_3);
            this.Controls.Add(this.launch_name_2);
            this.Controls.Add(this.launch_name_1);
            this.Controls.Add(this.script_file_label);
            this.Controls.Add(this.executable_select_button);
            this.Controls.Add(this.flexPDE_location);
            this.Controls.Add(this.flexpde_filename);
            this.Controls.Add(this.error_label);
            this.Controls.Add(this.s_label);
            this.Controls.Add(this.Time_label);
            this.Controls.Add(this.time_text);
            this.Controls.Add(this.timetostart_progressbar);
            this.Controls.Add(this.Stop_button);
            this.Controls.Add(this.Start_button);
            this.Name = "Form1";
            this.Text = "FlexPDE Startup Tool";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button Start_button;
        private System.Windows.Forms.Button Stop_button;
        private System.Windows.Forms.ProgressBar timetostart_progressbar;
        private System.Windows.Forms.TextBox time_text;
        private System.Windows.Forms.Label Time_label;
        private System.Windows.Forms.Label s_label;
        private System.Windows.Forms.Label error_label;
        private System.Windows.Forms.TextBox flexpde_filename;
        private System.Windows.Forms.Label flexPDE_location;
        private System.Windows.Forms.Button executable_select_button;
        private System.Windows.Forms.OpenFileDialog openFileDialog1;
        private System.Windows.Forms.Label script_file_label;
        private System.Windows.Forms.TextBox launch_name_1;
        private System.Windows.Forms.TextBox launch_name_2;
        private System.Windows.Forms.TextBox launch_name_3;
        private System.Windows.Forms.TextBox launch_name_4;
        private System.Windows.Forms.TextBox launch_name_5;
        private System.Windows.Forms.CheckBox launch_check_1;
        private System.Windows.Forms.CheckBox launch_check_2;
        private System.Windows.Forms.CheckBox launch_check_3;
        private System.Windows.Forms.CheckBox launch_check_4;
        private System.Windows.Forms.CheckBox launch_check_5;
        private System.Windows.Forms.Button launch_select_1;
        private System.Windows.Forms.Button launch_select_2;
        private System.Windows.Forms.Button launch_select_3;
        private System.Windows.Forms.Button launch_select_4;
        private System.Windows.Forms.Button launch_select_5;
        private System.Windows.Forms.CheckBox run_check;
        private System.Windows.Forms.CheckBox exit_check;
        private System.Windows.Forms.CheckBox minimized_check;
        private System.Windows.Forms.CheckBox silent_check;
        private System.Windows.Forms.Button launch_select_6;
        private System.Windows.Forms.CheckBox launch_check_6;
        private System.Windows.Forms.TextBox launch_name_6;
        private System.Windows.Forms.Button logout_button;
        private System.Windows.Forms.TextBox session_no_text;
        private System.Windows.Forms.Label session_id_label;
    }
}

