using System;
using System.Linq;
using System.Numerics;
using System.Collections;
using System.Data;
using System.IO;


namespace GAfitting
{
    public class TabulatedMaterial
    {
        /* ----------------------------------------------------------------- */
        /* |                             accessors                         | */
        /* ----------------------------------------------------------------- */

        //array storing the tabulated wavelengths
        public ArrayList Lambda { get; private set; }
        //array storing complex permittivities
        public ArrayList Epsilon { get; private set; }
        //units in which the wavelength is stored in the file
        public FDTDutils.units LambdaUnits { get; set; }
        //are we using N,K values or real, imaginary (eps) values
        public Boolean isNK { get; set; }
        //separator caracter used in source file
        public string Separator { get; set; }



        /* ----------------------------------------------------------------- */
        /* |                          public methods                       | */
        /* ----------------------------------------------------------------- */

        public TabulatedMaterial()
        {
            Lambda = new ArrayList();
            Epsilon = new ArrayList();
            LambdaUnits = FDTDutils.units.um;
        }

        //Load a tabulated lambda, n, k data file and fills the arrays
        public bool LoadDataFile(string filename)
        {
            try
            {
                //datatable is probably overkill here...     
                DataTable table = new DataTable();
                table.Columns.Add("lambda");
                table.Columns.Add("val1");
                table.Columns.Add("val2");

                var lines = File.ReadAllLines(filename).ToList();
                lines.ForEach(line => table.Rows.Add(line.Split(new String[] {Separator}, StringSplitOptions.RemoveEmptyEntries)));

                ConvertTable(table);
            }
            catch (Exception e)
            {
                Console.WriteLine("The file could not be read:");
                Console.WriteLine(e.Message);
                return false;
            }
            return true;
        }

        //returns the dielectric function at a specified lambda
        public virtual Complex GetEpsilon(double lambda, params double[] values)
        {
            int idx = Lambda.IndexOf(lambda);
            if (idx > -1)
                return (Complex)Epsilon[idx];
            else
                return (Complex)0.0;
        }

        //writes the material definition to an XML file 
        public virtual void WriteXML(string folder, string matname, double[] result)
        {
            // not implemented, must override 
        }


        /* ----------------------------------------------------------------- */
        /* |                       private members                         | */
        /* ----------------------------------------------------------------- */

        //Converts the datatable into epsilon and lambda arrays
        private void ConvertTable(DataTable table)
        {

            if (LambdaUnits == FDTDutils.units.Hz)
            {
                DataView dv = table.DefaultView;
                dv.Sort = "lambda DESC"; //we want ascending wavelengths, so we set it to descending frequencies
                table = dv.ToTable();
            }

            foreach (DataRow row in table.Rows)
            {
                double val1 = Convert.ToDouble(row["val1"].ToString());
                double val2 = Convert.ToDouble(row["val2"].ToString());
                double lambda = Convert.ToDouble(row["lambda"].ToString());

                Lambda.Add(lambda);
                Complex eps;
                if (isNK)
                    eps = ConvertNKToEpsilon(val1, val2);
                else
                    eps = ConvertRealImagToEpsilon(val1, val2);

                Epsilon.Add(eps);
            }

        }
        //converts the n,k values into a complex dielectric constant espilon
        private Complex ConvertNKToEpsilon(double n, double k)
        {
            return new Complex(n * n - k * k, 2 * n * k);
        }

        //converts the real,imaginary (eps) values into a complex dielectric 
        //constant espilon
        private Complex ConvertRealImagToEpsilon(double realeps, double imageps)
        {
            return new Complex(realeps, imageps);
        }
      

    }
}
