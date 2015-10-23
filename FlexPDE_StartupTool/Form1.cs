using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Diagnostics;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace FlexPDE_StartupTool
{
    public partial class Form1 : Form
    {
        const double time_default = 60000.0;
        double wait_time;
        bool counting_down = false;

        public Form1()
        {
            InitializeComponent();
            InitializeBackgroundWorker();
            Reset_Time();
        }

        // Set up the BackgroundWorker object by 
        // attaching event handlers. 
        private void InitializeBackgroundWorker()
        {
            backgroundworker.DoWork += new DoWorkEventHandler(Worker_Sleep);
            backgroundworker.RunWorkerCompleted += new RunWorkerCompletedEventHandler(Close_Worker);
            backgroundworker.ProgressChanged += new ProgressChangedEventHandler(Progress_Update);
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        BackgroundWorker backgroundworker;
        DateTime start_time;
        DateTime end_time;
        private void start_Click(object sender, EventArgs e)
        {
            if (!double.TryParse(time_text.Text, out wait_time))
            {
                MessageBox.Show("Error - Could not parse start time");
                return;
            }

            // convert initial time from seconds to milliseconds
            wait_time *= 1000.0;
            counting_down = true;
            error_label.Text = "";

            start_time = DateTime.Now;
            end_time = start_time + new TimeSpan(0, 0, 0, 0, (int)Math.Round(wait_time));
            backgroundworker.RunWorkerAsync(wait_time);
        }

        private void stop_Click(object sender, EventArgs e)
        {
            backgroundworker.CancelAsync();
        }

        private void Reset_Time()
        {
            wait_time = time_default;
            time_text.Text = (time_default / 1000.0).ToString("F1");
            time_text.Update();
            timetostart_progressbar.Value = 0;
        }

        void Worker_Sleep(object sender, DoWorkEventArgs e)
        {
            BackgroundWorker worker = sender as BackgroundWorker;

            e.Result = Worker_Wait((double)e.Argument, worker, e);

            if ((bool)e.Result == true)
                Launch_FlexPDE();
        }

        bool Worker_Wait(double time, BackgroundWorker worker, DoWorkEventArgs e)
        {
            while(counting_down)
            {
                if (worker.CancellationPending)
                {
                    e.Cancel = true;
                    counting_down = false;
                }

                Thread.Sleep(100);
                int progress = (int)Math.Round(100.0 * (DateTime.Now - start_time).TotalMilliseconds / time);
                if (progress > 100)
                    progress = 100;
                worker.ReportProgress(progress);
                if (progress > 99)
                    return true;
            }

            return false;
        }
        
        private void Progress_Update(object sender, ProgressChangedEventArgs e)
        {
            time_text.Text = ((end_time - DateTime.Now).TotalMilliseconds / 1000.0).ToString("F1");
            time_text.Update();
            timetostart_progressbar.Value = e.ProgressPercentage;
        }

        void Close_Worker(object sender, RunWorkerCompletedEventArgs e)
        {
            // First, handle the case where an exception was thrown.
            if (e.Error != null)
            {
                MessageBox.Show(e.Error.Message);
            }

            counting_down = false;
            Reset_Time();
        }

        private void Launch_FlexPDE()
        {
            string args = "";
            if (run_check.Checked)
                args += "-R ";
            if (exit_check.Checked)
                args += "-X ";
            if (minimized_check.Checked)
                args += "-M ";
            if (silent_check.Checked)
                args += "-S ";

            string orig_dir = Directory.GetCurrentDirectory();

            if (launch_check_1.Checked && launch_name_1.Text != "") 
                Start_FlexPDE_Process(args, orig_dir, launch_name_1.Text);
            if (launch_check_2.Checked && launch_name_2.Text != "")
                Start_FlexPDE_Process(args, orig_dir, launch_name_2.Text);
            if (launch_check_3.Checked && launch_name_3.Text != "")
                Start_FlexPDE_Process(args, orig_dir, launch_name_3.Text);
            if (launch_check_4.Checked && launch_name_4.Text != "")
                Start_FlexPDE_Process(args, orig_dir, launch_name_4.Text);
            if (launch_check_5.Checked && launch_name_5.Text != "")
                Start_FlexPDE_Process(args, orig_dir, launch_name_5.Text);
            if (launch_check_6.Checked && launch_name_6.Text != "")
                Start_FlexPDE_Process(args, orig_dir, launch_name_6.Text);
        }

        private void Start_FlexPDE_Process(string args, string orig_dir, string target_dir)
        {
            Directory.SetCurrentDirectory(Path.GetDirectoryName(target_dir));
            Process.Start(flexpde_filename.Text, args + Path.GetFileName(target_dir));
            Directory.SetCurrentDirectory(orig_dir);
            Thread.Sleep(1000);
        }

        private void select_Click(object sender, EventArgs e)
        {
            // Show the "Make new folder" button
            openFileDialog1.InitialDirectory = "C:\\FlexPDE6";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                flexpde_filename.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }

        private void launch_select_1_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = "C:\\Users\\Edmund\\Documents\\Visual Studio 2013\\Projects\\Iterative_Greens_Function\\TwoD_ThomasFermiPoisson\\bin\\x64\\";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                launch_name_1.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }
        private void launch_select_2_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = "C:\\Users\\Edmund\\Documents\\Visual Studio 2013\\Projects\\Iterative_Greens_Function\\TwoD_ThomasFermiPoisson\\bin\\x64\\";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                launch_name_2.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }
        private void launch_select_3_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = "C:\\Users\\Edmund\\Documents\\Visual Studio 2013\\Projects\\Iterative_Greens_Function\\TwoD_ThomasFermiPoisson\\bin\\x64\\";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                launch_name_3.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }
        private void launch_select_4_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = "C:\\Users\\Edmund\\Documents\\Visual Studio 2013\\Projects\\Iterative_Greens_Function\\TwoD_ThomasFermiPoisson\\bin\\x64\\";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                launch_name_4.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }
        private void launch_select_5_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = "C:\\Users\\Edmund\\Documents\\Visual Studio 2013\\Projects\\Iterative_Greens_Function\\TwoD_ThomasFermiPoisson\\bin\\x64\\";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                launch_name_5.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }
        private void launch_select_6_Click(object sender, EventArgs e)
        {
            openFileDialog1.InitialDirectory = "C:\\Users\\Edmund\\Documents\\Visual Studio 2013\\Projects\\Iterative_Greens_Function\\TwoD_ThomasFermiPoisson\\bin\\x64\\";

            if (openFileDialog1.ShowDialog() == DialogResult.OK)
            {
                launch_name_6.Text = openFileDialog1.FileName;
            }
            else
            {
                // Selection canceled
                MessageBox.Show("Operation canceled.");
            }
        }

        private void logout_button_Click(object sender, EventArgs e)
        {
            int i;
            if (int.TryParse(session_no_text.Text, out i))
            {
                ProcessStartInfo newProcessInfo = new ProcessStartInfo();
                newProcessInfo.FileName = @"C:\Windows\System32\tscon.exe";
                newProcessInfo.Verb = "runas";

                string cmd_text = i.ToString() + " /dest:console";
                newProcessInfo.Arguments = cmd_text;

                Process.Start(newProcessInfo);
            }
            else
                MessageBox.Show("Error - Invalid session ID");
        }
    }
}

