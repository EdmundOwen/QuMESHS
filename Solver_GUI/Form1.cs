using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Solver_GUI
{
    public partial class Solver_GUI : Form
    {
        int dimension = 2;

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
    }
}
