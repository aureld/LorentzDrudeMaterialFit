using System;
using System.IO;
using System.Numerics;
using System.Xml;

namespace GAfitting
{
    public class TabulatedLorentzDrude : TabulatedMaterial, IGeneticAlgorithm
    {

        /* ----------------------------------------------------------------- */
        /* |                             accessors                         | */
        /* ----------------------------------------------------------------- */

        public double[] UpperBounds { get; set; } //Upper constraints for the Lorentz-Drude coefficients
        public double[] LowerBounds { get; set; } //Lower constraints for the Lorentz-Drude coefficients
        public bool FitEpsInf { get; set; } //do we want to free the eps_inf coefficient? if no, eps_inf = 1.0

        /* ----------------------------------------------------------------- */
        /* |                          public methods                       | */
        /* ----------------------------------------------------------------- */

        public TabulatedLorentzDrude()
            : base()
        {
            FitEpsInf = false;
        }

        public TabulatedLorentzDrude(int nbterms, bool fitepsinf, double[] low, double[] up)
            : base()
        {
            FitEpsInf = fitepsinf;

            double[] lbounds, ubounds;

            if (fitepsinf)
            {
                lbounds = new double[nbterms * 3 + 1];
                ubounds = new double[nbterms * 3 + 1];
                //if we want to fit eps_inf, we add a term in front
                lbounds[0] = low[0];
                ubounds[0] = up[0]; 
                for (int i = 1; i < nbterms * 3 + 1; i += 3)
                {
                    lbounds[i] = low[1];
                    lbounds[i + 1] = low[2];
                    lbounds[i + 2] = low[3];
                    ubounds[i] = up[1];
                    ubounds[i + 1] = up[2];
                    ubounds[i + 2] = up[3];
                }
            }
            else
            {
                lbounds = new double[nbterms * 3];
                ubounds = new double[nbterms * 3];
                for (int i = 0; i < nbterms * 3; i += 3)
                {
                    lbounds[i] = low[0];
                    lbounds[i + 1] = low[1];
                    lbounds[i + 2] = low[2];
                    ubounds[i] = up[0];
                    ubounds[i + 1] = up[1];
                    ubounds[i + 2] = up[2];
                }
            }

            LowerBounds = lbounds;
            UpperBounds = ubounds;

        }

        //Fitness score of the current candidates dataset against the tabulated data
        public double FitnessFunction(double[] candidates)
        {
            double fit = 0.0;
            int idx = 0;

            //scale candidate genes
            double[] scaledGenes = ScaleGenes(candidates);

            foreach (double lambda in Lambda)
            {
                //calculate the dielectric constant based on candidate dataset
                Complex calculatedEpsilon = GetEpsilon(lambda, scaledGenes);
                //compare the calculated and measured values 
                Complex measuredEpsilon = (Complex)Epsilon[idx];
                Complex delta = calculatedEpsilon - measuredEpsilon;
                fit += delta.Magnitude; //we are looking at minimizing delta
                idx++;
            }
            return 1 / (1 + fit); //we want to maximize fitness
        }

        //returns the estimated value of the dielectric constant for the current set of genes
        public override Complex GetEpsilon(double lambda, params double[] values)
        {
            if (values.GetLength(0) < 3)
                throw new ArgumentOutOfRangeException("Not enough arguments.");

            double omega;
            //convert the wavelength to angular frequency
            if (LambdaUnits != FDTDutils.units.Hz)
                omega = FDTDutils.LambdaToOmega(lambda, LambdaUnits);
            else
                omega = lambda;

            Complex eps = 1.0; // default to air

            //overrides default to allow eps_inf as a coefficient
            if (FitEpsInf)
                eps = values[0]; // epsilon_inf

            // grouping of values:
            //values[1] = S1, values[2] = omega_1, values[3] = gamma_1 
            //values[4] = S2, values[5] = omega_2, values[6] = gamma_2  etc.

            //special case for Drude term
            double s_i = values[0];
            double omega_p = values[1];
            double gamma_i = values[2];

            if (FitEpsInf)
            {
                s_i = values[1];
                omega_p = values[2];
                gamma_i = values[3];
            }

            double omega_i = 0.0; // no resonance for Drude term

            eps += -s_i * omega_p * omega_p / (Complex.ImaginaryOne * omega * gamma_i + omega * omega);

            //Lorentz terms
            int imin = 3;
            if (FitEpsInf)
                imin = 4;

            for (int i = imin; i < values.GetLength(0); i += 3)
            {
                s_i = values[i];
                omega_i = values[i + 1];
                gamma_i = values[i + 2];
                eps += s_i * omega_p * omega_p / (omega_i * omega_i - omega * omega - Complex.ImaginaryOne * omega * gamma_i);
            }

            return eps;
        }

        //returns the average error of the fitting -- TODO: fix mismatch matlab / this 
        public double GetFittingError(params double[] values)
        {
            double error = 0.0;
            foreach (double lambda in Lambda)
            {
                //calculated epsilon
                double epsRcalc = GetEpsilon(lambda, values).Real;
                double epsIcalc = GetEpsilon(lambda, values).Imaginary;

                //reference epsilon
                double epsR = base.GetEpsilon(lambda).Real;
                double epsI = base.GetEpsilon(lambda).Imaginary;

                //Error calculation
                
                double epsRnormErr = Math.Pow(Math.Sqrt(Math.Abs(epsRcalc - epsR)), 2);
                double epsInormErr = Math.Pow(Math.Sqrt(Math.Abs(epsIcalc - epsI)), 2);
                double epsRnorm = Math.Abs(epsR);
                double epsInorm = Math.Abs(epsI);

                error += (epsRnormErr + epsInormErr) / (epsRnorm + epsInorm);
                
            }

            return error;
        }

        //scale candidate genes using upper and lower bounds 
        public double[] ScaleGenes(double[] candidates)
        {
            double[] result = new double[candidates.GetLength(0)];
            for (int i = 0; i < result.GetLength(0); i++)
            {
                double range = UpperBounds[i] - LowerBounds[i];
                result[i] = LowerBounds[i] + candidates[i] * range;
            }
            return result;
        }

        //writes the material definition to an XML file 
        public override void WriteXML(string folder, string matname, double[] result)
        {
            using (XmlWriter xw = XmlWriter.Create(Path.Combine(folder, matname + ".xml")))
            {
                xw.WriteStartDocument();
                xw.WriteStartElement("ProfileDesignerData");
                xw.WriteAttributeString("optiwave_xml_file_version", "1.0");
                xw.WriteStartElement("MatDielectricLorentzDrude");
                xw.WriteAttributeString("name", matname);
                if (FitEpsInf)
                {
                    xw.WriteStartElement("Einf");
                    xw.WriteAttributeString("X", result[0].ToString());
                    xw.WriteAttributeString("Y", result[0].ToString());
                    xw.WriteAttributeString("Z", result[0].ToString());
                    xw.WriteEndElement();
                    //Drude coefs
                    xw.WriteStartElement("LDTerm");
                    xw.WriteAttributeString("S", result[1].ToString());
                    xw.WriteAttributeString("P", result[2].ToString());
                    xw.WriteAttributeString("R", "0.0");
                    xw.WriteAttributeString("D", result[3].ToString());
                    xw.WriteEndElement();
                    //Lorentz coefs
                    for (int i = 4; i < result.Length; i += 3)
                    {
                        xw.WriteStartElement("LDTerm");
                        xw.WriteAttributeString("S", result[i].ToString());
                        xw.WriteAttributeString("P", result[2].ToString());
                        xw.WriteAttributeString("R", result[i + 1].ToString());
                        xw.WriteAttributeString("D", result[i + 2].ToString());
                        xw.WriteEndElement();
                    }
                }
                else
                {
                    xw.WriteStartElement("Einf");
                    xw.WriteAttributeString("X", "1.0");
                    xw.WriteAttributeString("Y", "1.0");
                    xw.WriteAttributeString("Z", "1.0");
                    xw.WriteEndElement();
                    //Drude coefs
                    xw.WriteStartElement("LDTerm");
                    xw.WriteAttributeString("S", result[0].ToString());
                    xw.WriteAttributeString("P", result[1].ToString());
                    xw.WriteAttributeString("R", "0.0");
                    xw.WriteAttributeString("D", result[2].ToString());
                    xw.WriteEndElement();
                    //Lorentz coefs
                    for (int i = 3; i < result.Length; i += 3)
                    {
                        xw.WriteStartElement("LDTerm");
                        xw.WriteAttributeString("S", result[i].ToString());
                        xw.WriteAttributeString("P", result[1].ToString());
                        xw.WriteAttributeString("R", result[i + 1].ToString());
                        xw.WriteAttributeString("D", result[i + 2].ToString());
                        xw.WriteEndElement();
                    }
                }

                xw.WriteEndElement();
                xw.WriteEndElement();
                xw.WriteEndDocument();
            }
        }


        /* ----------------------------------------------------------------- */
        /* |                       private members                         | */
        /* ----------------------------------------------------------------- */


       
    }
}
