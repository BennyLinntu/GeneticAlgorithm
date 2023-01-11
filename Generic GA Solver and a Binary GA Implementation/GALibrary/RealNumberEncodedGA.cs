using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GALibrary
{
    // data fields
    public enum RealNumberCrossoverType
    {
        LVD, SVD, MOES, TES, FBMS
    }
    public enum RealNUmberMutationType
    {
        DynamicMutation
    }


    class RealNumberEncodedGA : GenericSolver<double>
    {

        double[] lowerBound;
        double[] upperBound;
        double degreeofNonnuiFormity = 1;
        public RealNumberCrossoverType crossoverType { set; get; } = RealNumberCrossoverType.LVD;
        public RealNUmberMutationType mutationType { set; get; } = RealNUmberMutationType.DynamicMutation;



        public RealNumberEncodedGA(int numberofGenes, double[] lowerBound, double[] upperBound, GAOptimizationType optimizationType, ObjectiveFunction<double> objectiveFunction) : base(numberofGenes, optimizationType, objectiveFunction)
        {
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }
        public override void InitializePopulation()
        {
            for (int row = 0; row < PopulationSize; row++)
            {
                for (int column = 0; column < numberofGenes; column++)
                {
                    Chromosomes[row][column] = lowerBound[column] + rnd.NextDouble() * (upperBound[column] - lowerBound[column]);
                }
                objectiveValue[row] = objectivefunction(Chromosomes[row]);
            }

        }

        



    }
}
