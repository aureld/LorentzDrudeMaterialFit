// tester for the genetic algorithm fitting 

using System;
using GAfitting;
using System.IO;
using System.Xml;
using System.Diagnostics;



public class Test
{
	public static void Main()
	{
        double eh = FDTDutils.ehbar;

        /* ----------------------------------------------------------------- */
        /* |                          PARAMETERS                           | */
        /* ----------------------------------------------------------------- */

        //boundaries of the search space (strength, omega (rad/s), damping (rad/s))
        //add nmin nmax in front if fitepsinf = true
        double[] up = { 10.0, 100 * eh, 100 * eh };
        double[] low = { 0.001, 0.001 * eh, 0.001 * eh };
        Boolean fitepsinf = false; // do we want to free the eps_inf coefficient? if no, eps_inf = 1.0
        int nb_terms = 6; //Number of LD terms we want (including the Drude term)
        double co = 0.8;            //  Crossover rate
        double mut = 0.05;          //  Mutation rate
        int pop_size = 100;         //  Population size
        int generations = 50000;    //  Number of generations
        Boolean elitism = true;     //do we want to keep the best of each generation?
        Boolean nk = true; //is the file in (n,k) form? if not, (real(eps),imag(eps)) is assumed 
        string separator = " "; //column separator
        FDTDutils.units units = FDTDutils.units.um; //wavelength units of the input file 
        string folder = @".";
        string matname = "Silver_Johnson";
        string extension = ".txt";

        /* ----------------------------------------------------------------- */
        /* |                          main program                         | */
        /* ----------------------------------------------------------------- */

        
        TabulatedLorentzDrude diel = new TabulatedLorentzDrude(nb_terms, fitepsinf, low, up);
        diel.LambdaUnits =units; 
        diel.isNK = nk; 
        diel.Separator = separator;
        diel.LoadDataFile(Path.Combine(folder,matname + extension));

        int genome = fitepsinf? 3*nb_terms+1 : 3*nb_terms; //  Genome size
        GA ga = new GA(co, mut, pop_size, generations, genome);
		
        //GA algo init
        ga.FitnessFunction = new GAFunction(diel.FitnessFunction);
        ga.FitnessFile = Path.Combine(folder,matname + "_fitness.csv"); //fitness score output to file
        ga.Elitism = elitism;

        Stopwatch stopwatch = new Stopwatch();
        stopwatch.Start();
        ga.Go(); 
        stopwatch.Stop();
        Console.WriteLine("elapsed time (ms): {0}", stopwatch.ElapsedMilliseconds);

        //returns best fitness score and relative error <-- TODO fix relative error vs matlab
        double[] values;
        double fitness;
        ga.GetBest(out values, out fitness);
        double[] result = diel.ScaleGenes(values);
        System.Console.WriteLine("Best ({0}), error {1}", fitness, diel.GetFittingError(result));

        //save coefs to text
        using (StreamWriter sw = new StreamWriter(Path.Combine(folder,matname + "_coefs.txt")))
        {
            for (int i = 0; i < result.Length; i++)
                sw.WriteLine("{0} ", result[i]);
        }

        //save material LD coefs to XML
        diel.WriteXML(folder, matname, result);

    }
}
