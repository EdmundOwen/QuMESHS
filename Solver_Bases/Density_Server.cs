using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Solver_Bases.Layers;
using NetworkCommsDotNet;

namespace Solver_Bases
{
    /// <summary>
    /// A class which acts as the server side of an externally calculated density
    /// </summary>
    public class Density_Server : Density_Base
    {
        string clientIP;
        int clientPort;

        SpinResolved_Data tmp_density;

        public Density_Server(double temperature, double dz, int nz, double zmin)
            : base(temperature, 1.0, 1.0, dz, 1, 1, nz, 0.0, 0.0, zmin)
        {
            Create_Connection();
        }

        public Density_Server(double temperature, double dy, double dz, int ny, int nz, double ymin, double zmin)
            : base(temperature, 1.0, dy, dz, 1, ny, nz, 0.0, ymin, zmin)
        {
            Create_Connection();
        }

        public Density_Server(double temperature, double dx, double dy, double dz, int nx, int ny, int nz, double xmin, double ymin, double zmin) 
            : base(temperature, dx, dy, dz, nx, ny, nz, xmin, ymin, zmin)
        {
            Create_Connection();
        }

        public override void Get_ChargeDensity(ILayer[] layers, ref SpinResolved_Data density, Band_Data chem_pot)
        {
            // send data to client
            NetworkComms.SendObject("layers", clientIP, clientPort, layers);
            NetworkComms.SendObject("density", clientIP, clientPort, density);
            NetworkComms.SendObject("chem_pot", clientIP, clientPort, chem_pot);

            // Call the save_density routine
            //NetworkComms.AppendGlobalIncomingPacketHandler<SpinResolved_Data>("new_density", Save_Density);

            //Start listening for incoming connections
            TCPConnection.StartListening(true);

            //Print out the IPs and ports we are now listening on
            Console.WriteLine("Server listening for TCP connection on:");
            foreach (System.Net.IPEndPoint localEndPoint in TCPConnection.ExistingLocalListenEndPoints()) Console.WriteLine("{0}:{1}", localEndPoint.Address, localEndPoint.Port);
 
            throw new NotImplementedException();
        }

        public override double Get_Chemical_Potential(double x, double y, double z, ILayer[] layers, double temperature_input)
        {
            throw new NotImplementedException();
        }

        private void Save_Density()
        {
        }

        private void Create_Connection()
        {
            //Request client IP and port number
            Console.WriteLine("Please enter the client IP and port in the format 192.168.0.1:10000 and press return:");
            string clientInfo = Console.ReadLine();
            
            //Parse the necessary information out of the provided string
            clientIP = clientInfo.Split(':').First();
            clientPort = int.Parse(clientInfo.Split(':').Last());

 
            throw new NotImplementedException();
        }

        public override void Close()
        {
            NetworkComms.Shutdown();
            Console.WriteLine("Closing density solver");
        }
    }
}
