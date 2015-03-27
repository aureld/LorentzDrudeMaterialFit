
namespace GAfitting
{
    interface IGeneticAlgorithm
    {
        double FitnessFunction(double[] candidates);
        double GetFittingError(params double[] values);
    }
}
