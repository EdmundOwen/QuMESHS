using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Geometry;

namespace Solver_Bases.Layers
{
    public interface ILayer
    {
        bool Dopents_Frozen_Out(double temperature);
        void Set_Dopents(double acceptor_concentration, double donor_concentration);

        #region Access_Geometry
        bool InLayer(double z);
        bool InLayer(double y, double z);
        bool InLayer(double x, double y, double z);
        Geometry_Type Geometry { get; }
        ILayer Get_Component(int component_no);
        int No_Components { get; }

        double Xmin { get; }
        double Xmax { get; }
        double Ymin { get; }
        double Ymax { get; }
        double Zmin { get; }
        double Zmax { get; }
        #endregion

        #region Layer_Properties
        int Layer_No { get; }
        Material Material { get; }
        double Permitivity { get; }
        double Band_Gap { get; }
        double Donor_Energy { get; }
        double Acceptor_Energy { get; }
        double Donor_Conc { get; }
        double Acceptor_Conc { get; }
        double Dopent_FreezeOut_T { get; }
        #endregion
    }
}
