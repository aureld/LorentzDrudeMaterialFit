using System;

namespace GAfitting
{
    public static class FDTDutils
    {
        //Converts a wavelength in um to an angular frequency in rad/s
        public static double LambdaToOmega(double lambda, FDTDutils.units units)
        {
            double factor = 0.0;
            switch (units)
            {
                case FDTDutils.units.nm:
                    factor = 1e9;
                    break;
                case FDTDutils.units.um:
                    factor = 1e6;
                    break;
                default:
                    factor = 1.0;
                    break;
            }


            return (2 * Math.PI * c_light * factor / lambda);
        }

        public static double c_light = 299792458.0; //speed of light in vacuum 
        public static double ehbar = 1.51926751447914e+015; // e / hbar with hbar = h/2pi and e=1.6e-19 


        public enum units { um, nm, Hz }; // wavelength / frequency units


    }
}
