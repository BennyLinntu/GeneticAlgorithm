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

        private void LVD_Crossover(int father, int mother, int child_a, int child_b)
        {
            for (int i = 0; i < numberofGenes; i++)
            {
                double b_large;
                double alpha = rnd.NextDouble();
                if (chromosomes[father][i] > chromosomes[mother][i])
                    b_large = chromosomes[father][i];
                else
                    b_large = chromosomes[mother][i];

                chromosomes[child_a][i] = alpha * lowerBound[i] + (1 - alpha) * b_large;
                chromosomes[child_b][i] = alpha * b_large + (1 - alpha) * upperBound[i];
            }

        }
        private void SVD_Crossover(int father, int mother, int child_a, int child_b)
        {
            for (int i = 0; i < numberofGenes; i++)
            {
                double b_small;
                double alpha = rnd.NextDouble();
                if (chromosomes[father][i] < chromosomes[mother][i])
                    b_small = chromosomes[father][i];
                else
                    b_small = chromosomes[mother][i];

                chromosomes[child_a][i] = alpha * lowerBound[i] + (1.0 - alpha) * b_small;
                chromosomes[child_b][i] = alpha * b_small + (1.0 - alpha) * upperBound[i];
            }

        }
        private void MOES_Crossover(int father, int mother, int child_a, int child_b)
        {
            for (int i = 0; i < numberofGenes; i++)
            {
                double b_small;
                double b_large;
                double alpha = rnd.NextDouble();
                double u = rnd.NextDouble();
                if (chromosomes[father][i] < chromosomes[mother][i])
                {
                    b_small = chromosomes[father][i];
                    b_large = chromosomes[mother][i];
                }
                else
                {
                    b_small = chromosomes[mother][i];
                    b_large = chromosomes[father][i];
                }

                chromosomes[child_a][i] = alpha * b_small + (1.0 - alpha) * b_large;
                if (u > 0.5)
                    chromosomes[child_b][i] = alpha * b_large + (1.0 - alpha) * upperBound[i];
                else
                    chromosomes[child_b][i] = alpha * lowerBound[i] + (1.0 - alpha) * b_small;
            }

        }

        private void TES_Crossover(int father, int mother, int child_a, int child_b)
        {
            for (int i = 0; i < numberofGenes; i++)
            {
                double b_small;
                double b_large;
                double alpha = rnd.NextDouble();
                double u = rnd.NextDouble();
                if (chromosomes[father][i] < chromosomes[mother][i])
                {
                    b_small = chromosomes[father][i];
                    b_large = chromosomes[mother][i];
                }
                else
                {
                    b_small = chromosomes[mother][i];
                    b_large = chromosomes[father][i];
                }

                chromosomes[child_a][i] = alpha * lowerBound[i] + (1.0 - alpha) * b_small;
                chromosomes[child_b][i] = alpha * b_large + (1.0 - alpha) * upperBound[i];
            }

        }

        private void FBMS_Crossover(int father, int mother, int child_a, int child_b)
        {
            for (int i = 0; i < numberofGenes; i++)
            {
                double b_small;
                double b_large;
                double alpha = rnd.NextDouble();
                if (chromosomes[father][i] < chromosomes[mother][i])
                {
                    b_small = chromosomes[father][i];
                    b_large = chromosomes[mother][i];
                }
                else
                {
                    b_small = chromosomes[mother][i];
                    b_large = chromosomes[father][i];
                }

                chromosomes[child_a][i] = alpha * b_small + (1.0 - alpha) * b_large;
                chromosomes[child_b][i] = alpha * b_large + (1.0 - alpha) * b_small;
            }

        }



    }
}
