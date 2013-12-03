using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public abstract class Layer : ILayer
    {
        // layer information
        IGeom geom;
        protected Material material;
        protected double permitivity;
        int layer_no;

        // material parameters
        protected double band_gap;
        protected bool allow_donors = false;        // by default, don't allow donors to be set
        protected double donor_energy;
        protected double acceptor_energy;
        protected double donor_concentration = 0.0;
        protected double acceptor_concentration = 0.0;
        protected double freeze_out_temp = 0.0;     // default for Layer is no freezing out of donors
                
        public Layer(IGeom Geom, int layer_no)
        {
            this.layer_no=layer_no;
            this.geom = Geom;
            Set_Material_Parameters();
        }

        protected abstract void Set_Material_Parameters();
        public void Set_Dopents(double acceptor_concentration, double donor_concentration)
        {
            if (!allow_donors)
                throw new NotImplementedException("Error - you are not allowed to add donors to layer no " + Layer_No.ToString());
            this.acceptor_concentration = acceptor_concentration; this.donor_concentration = donor_concentration;
            Set_Freeze_Out_Temperature();
        }
        public double Permitivity { get { return permitivity; } }
        public Material Material { get { return material; } }
        internal abstract void Set_Freeze_Out_Temperature();
        internal void Set_Freeze_Out_Temperature(double freeze_out_T)
        {
            if (!allow_donors)
                throw new NotImplementedException("Error - you are not allowed to add donors to layer no " + Layer_No.ToString());
            this.freeze_out_temp = freeze_out_T;
        }

        /// <summary>
        /// returns true if the input temperature is less than the temperature at which the donors are frozen out
        /// </summary>
        public bool Dopents_Frozen_Out(double temperature)
        {
            return temperature < freeze_out_temp;
        }

        #region Access_Geometry
        public bool InLayer(double z)
        { return geom.InLayer(z); }
        public bool InLayer(double y, double z)
        { return geom.InLayer(y, z); }
        public bool InLayer(double x, double y, double z)
        { return geom.InLayer(x, y, z); }

        public double Xmin { get { return geom.Xmin; } }
        public double Xmax { get { return geom.Xmax; } }
        public double Ymin { get { return geom.Ymin; } }
        public double Ymax { get { return geom.Ymax; } }
        public double Zmin { get { return geom.Zmin; } }
        public double Zmax { get { return geom.Zmax; } }
        #endregion

        #region Layer_Properties
        public int Layer_No { get { return layer_no; } }
        public double Band_Gap { get { return band_gap; } }
        public double Donor_Energy { get { return donor_energy; } }
        public double Acceptor_Energy { get { return acceptor_energy; } }
        public double Donor_Conc { get { return donor_concentration; } }
        public double Acceptor_Conc { get { return acceptor_concentration; } }
        public double Dopent_FreezeOut_T { get { return freeze_out_temp; } }
        #endregion
    }
}