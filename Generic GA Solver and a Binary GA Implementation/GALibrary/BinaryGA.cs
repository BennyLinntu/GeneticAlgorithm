using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GALibrary
{
    public enum BinaryCrossoverType
    {
        OnePointCut, TwoPointCut, NPointCut
    }

    public class BinaryGA : GenericSolver<byte>
    {
        // data field
        int[] cutPoints;
        BinaryCrossoverType crossoverType = BinaryCrossoverType.OnePointCut;
        // Property
        public BinaryCrossoverType CrossoverType { get => crossoverType; set => crossoverType = value; }

        public BinaryGA(int numberofGenes, GAOptimizationType optimizationType, ObjectiveFunction<byte> objectiveFunction, BinaryCrossoverType crossoverType) : base(numberofGenes, optimizationType, objectiveFunction)
        {
            this.crossoverType = crossoverType;
            cutPoints = new int[numberofGenes];
        }

        // function

        double GetByteValueasDouble(byte b)
        {
            if(b == 1)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
        bool ReturnReversedFlag(bool flag)
        {
            if(flag)
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        public override void InitializePopulation()
        {
            for (int row = 0; row < PopulationSize; row++)
            {
                for (int column = 0; column < numberofGenes; column++)
                {
                    Chromosomes[row][column] = (byte)rnd.Next(2);
                }
                objectiveValue[row] = objectivefunction(Chromosomes[row]);
            }
        }
        public override int[] ReturnChromosomesViolation(byte[] chrom)
        {
            int numberofJobs = Convert.ToInt32(Math.Sqrt(numberofGenes));
            int[] violations = new int[2 * numberofJobs];
            int counts;
            for (int r = 0; r < numberofJobs; r++)
            {
                counts = 0;
                for (int c = 0; c < numberofJobs; c++)
                {
                    counts += chrom[r * numberofJobs + c];
                }
                violations[r] = Math.Abs(counts - 1);
            }

            for (int c = 0; c < numberofJobs; c++)
            {
                counts = 0;
                for (int r = 0; r < numberofJobs; r++)
                {
                    counts += chrom[r * numberofJobs + c];     
                }
                violations[numberofJobs + c] = Math.Abs(counts - 1);
            }

            return violations;
        }
        public override void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            for (int i = 0; i < numberofGenes; i++)
            {
                Chromosomes[aftermutation][i] = Chromosomes[beforemutation][i];
                if(mutatedFlag[i])
                {
                    if(GetByteValueasDouble(Chromosomes[beforemutation][i]) == 1)
                    {
                        Chromosomes[aftermutation][i] = 0;
                    }
                    else
                    {
                        Chromosomes[aftermutation][i] = 1;
                    }
                }
            }
            objectiveValue[aftermutation] = objectivefunction(Chromosomes[aftermutation]);
        }


        public override void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb)
        {
            switch(crossoverType)
            {
                case BinaryCrossoverType.OnePointCut:
                    cutPoints[0] = rnd.Next(numberofGenes);
                    for (int j = 0; j < numberofGenes; j++)
                    {
                        if(j < cutPoints[0])
                        {
                            chromosomes[childa][j] = chromosomes[father][j];
                            chromosomes[childb][j] = chromosomes[mother][j];
                        }
                        else
                        {
                            chromosomes[childb][j] = chromosomes[father][j];
                            chromosomes[childa][j] = chromosomes[mother][j];
                        }
                    }
                    break;
                case BinaryCrossoverType.TwoPointCut:

                    cutPoints[0] = rnd.Next(numberofGenes);
                    cutPoints[1] = rnd.Next(numberofGenes);
                    Array.Sort(cutPoints, 0, 2);
                    for (int j = 0; j < numberofGenes; j++)
                    {
                        if (cutPoints[0] <= j && j < cutPoints[1])
                        {
                            Chromosomes[childa][j] = Chromosomes[father][j];
                            Chromosomes[childb][j] = Chromosomes[mother][j];
                        }
                        else
                        {
                            Chromosomes[childb][j] = Chromosomes[father][j];
                            Chromosomes[childa][j] = Chromosomes[mother][j];
                        }
                    }
                    break;
                case BinaryCrossoverType.NPointCut:
                    int numberofCuts = rnd.Next(numberofGenes);
                    bool flag = false;
                    cutPoints = Enumerable.Range(0, numberofGenes - 1).OrderBy(x => rnd.Next()).Take(numberofCuts).ToList().ToArray();

                    Array.Sort(cutPoints, 0, numberofCuts);
                    for (int j = 0; j < numberofGenes; j++)
                    {
                        if(cutPoints.Contains(j))
                        {
                            flag = ReturnReversedFlag(flag);
                        }
                        if(flag)
                        {
                            Chromosomes[childa][j] = Chromosomes[father][j];
                            Chromosomes[childb][j] = Chromosomes[mother][j];
                        }
                        else
                        {
                            Chromosomes[childb][j] = Chromosomes[father][j];
                            Chromosomes[childa][j] = Chromosomes[mother][j];
                        }
                    }
                    break;
                default:
                    break;
                    
            }
            objectiveValue[childa] = objectivefunction(chromosomes[childa]);
            objectiveValue[childb] = objectivefunction(chromosomes[childb]);
        }
        public override void ClearChromsomesandObjectiveArray()
        {
            for (int i = 0; i < 3 * PopulationSize; i++)
            {
                for (int j = 0; j < numberofGenes; j++)
                {
                    Chromosomes[i][j] = 0;
                }
                objectiveValue[i] = 0;
            }
        }
        

    }
}
