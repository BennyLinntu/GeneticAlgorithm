using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GALibrary
{
    // enum 
    public enum PermutationCrossoverType
    {
        ParialMappedX, OrderX, PositionBasedX, OrderBasedX, CycleX, SubtourExchange
    }
    public enum PermuatationMutationType
    {
        Inversion, Insertion, Displacement, ReciprocalExchange, Heuristic
    }

    public class PermutationGA : GenericSolver<int>
    {
        public PermutationCrossoverType CrossoverType { set; get; } = PermutationCrossoverType.ParialMappedX;
        public PermuatationMutationType PermuatationMutationType { set; get; } = PermuatationMutationType.Inversion;


        public PermutationGA(int numberofGenes, GAOptimizationType optimizationType, ObjectiveFunction<int> objectiveFunction, PermutationCrossoverType crossoverType, PermuatationMutationType mutationType) : base(numberofGenes, optimizationType, objectiveFunction)
        {
            this.CrossoverType = crossoverType;
            this.PermuatationMutationType = mutationType;
        }

        public int[] ShuffleaArray(int n, int take)
        {
            int[] arr = Enumerable.Range(0, n).OrderBy(x => rnd.Next()).Take(take).ToList().ToArray();
            return arr;
        }


        public override void InitializePopulation()
        {
            for (int i = 0; i < PopulationSize; i++)
            {
                for (int j = 0; j < numberofGenes; j++)
                {
                    Chromosomes[i] = ShuffleaArray(numberofGenes, numberofGenes);
                }
                ObjectiveValue[i] = objectivefunction(Chromosomes[i]);
            }
        }

        void ParialMappedCrossover(int father, int mother, int childa, int childb)
        {
            int[] twoPoints = ShuffleaArray(numberofGenes, 2);
            Array.Sort(twoPoints);
            int Pmin = twoPoints[0];
            int Pmax = twoPoints[1];
            int[] Dad = Chromosomes[father];
            int[] Mom = Chromosomes[mother];
            int swappingLength = Math.Abs(Pmax - Pmin);
            int[] Map = new int[numberofGenes];
            for (int i = 0; i < numberofGenes; i++)
            {
                Map[i] = -1;
            }

            for (int i = 0; i < Pmax; i++)
            {
                if (Dad[i] == Mom[i])
                {
                    continue;
                }
                if (Map[Dad[i]] == -1 && Map[Mom[i]] == -1)
                {
                    Map[Dad[i]] = Mom[i];
                    Map[Mom[i]] = Dad[i];
                }
                else if (Map[Dad[i]] == -1)
                {
                    Map[Dad[i]] = Map[Mom[i]];
                    Map[Map[Mom[i]]] = Dad[i];
                    Map[Mom[i]] = -2;
                }
                else if (Map[Mom[i]] == -1)
                {
                    Map[Mom[i]] = Map[Dad[i]];
                    Map[Map[Dad[i]]] = Mom[i];
                    Map[Dad[i]] = -2;
                }
                else
                {
                    Map[Map[Mom[i]]] = Map[Dad[i]];
                    Map[Map[Dad[i]]] = Map[Mom[i]];
                    Map[Dad[i]] = -3;
                    Map[Mom[i]] = -3;
                }
            }
            for (int i = 0; i < numberofGenes; i++)
            {
                if (i <= Pmax && i >= Pmin)
                {
                    Chromosomes[childb][i] = Dad[i];
                    Chromosomes[childa][i] = Mom[i];
                }
                else
                {
                    if (Map[Dad[i]] < 0)
                    {
                        Chromosomes[childa][i] = Dad[i];
                    }
                    else
                    {
                        Chromosomes[childa][i] = Map[Dad[i]];
                    }
                    if (Map[Mom[i]] < 0)
                    {
                        Chromosomes[childb][i] = Mom[i];
                    }
                    else
                    {
                        Chromosomes[childb][i] = Map[Mom[i]];
                    }
                }
            }
        }

        void OrderCrossover(int father, int mother, int childa, int childb)
        {
            int[] twoPoints = ShuffleaArray(numberofGenes, 2);
            Array.Sort(twoPoints);
            int Pmin = twoPoints[0];
            int Pmax = twoPoints[1];
            int swappingLength = Math.Abs(Pmax - Pmin);

            int[] swappingPartFromF = new int[swappingLength];
            int[] swappingPartFromM = new int[swappingLength];
            int[] remainpartforA = new int[numberofGenes - swappingLength];
            int[] remainpartforB = new int[numberofGenes - swappingLength];

            int r = 0;
            int rf = 0;
            int rm = 0;
            for (int i = 0; i < swappingLength; i++)
            {
                swappingPartFromF[i] = Chromosomes[father][Pmin + i];
                swappingPartFromM[i] = Chromosomes[mother][Pmin + i];
            }
            for (int i = 0; i < numberofGenes; i++)
            {
                if (swappingPartFromF.Contains(Chromosomes[mother][i]) == false)
                {
                    remainpartforA[rf] = Chromosomes[mother][i];
                    rf++;
                }
                if (swappingPartFromM.Contains(Chromosomes[father][i]) == false)
                {
                    remainpartforB[rm] = Chromosomes[father][i];
                    rm++;
                }
            }
            for (int i = 0; i < numberofGenes; i++)
            {
                if (i < Pmin || i >= Pmax)
                {
                    Chromosomes[childa][i] = remainpartforA[r];
                    Chromosomes[childb][i] = remainpartforB[r];
                }
                else
                {
                    Chromosomes[childa][i] = Chromosomes[father][i];
                    Chromosomes[childb][i] = Chromosomes[mother][i];

                }
            }

        }

        void PostitionBasedCrossover(int father, int mother, int childa, int childb)
        {
            int[] cutPoints = ShuffleaArray(numberofGenes, rnd.Next(1, numberofGenes));
            int swappingLength = cutPoints.Length;
            int[] swappingPartFromF = new int[swappingLength];
            int[] swappingPartFromM = new int[swappingLength];
            int[] remainpartforA = new int[numberofGenes - swappingLength];
            int[] remainpartforB = new int[numberofGenes - swappingLength];

            int s = 0;
            int r = 0;
            int rf = 0;
            int rm = 0;
            for (int i = 0; i < swappingLength; i++)
            {
                swappingPartFromF[i] = Chromosomes[father][cutPoints[i]];
                swappingPartFromM[i] = Chromosomes[mother][cutPoints[i]];
            }
            for (int i = 0; i < numberofGenes; i++)
            {
                if (swappingPartFromF.Contains(Chromosomes[mother][i]) == false)
                {
                    remainpartforA[rf] = Chromosomes[mother][i];
                    rf++;
                }
                if (swappingPartFromM.Contains(Chromosomes[father][i]) == false)
                {
                    remainpartforB[rm] = Chromosomes[father][i];
                    rm++;
                }
            }
            for (int i = 0; i < numberofGenes; i++)
            {

                if (cutPoints.Contains(i) == true)
                {
                    Chromosomes[childa][i] = swappingPartFromF[s];
                    Chromosomes[childb][i] = swappingPartFromM[s];
                    s++;
                }
                else
                {
                    Chromosomes[childa][i] = remainpartforA[r];
                    Chromosomes[childb][i] = remainpartforB[r];
                }

            }
        }

        void OrderBasedCrossover(int father, int mother, int childa, int childb)
        {
            int[] cutPoint_1 = ShuffleaArray(numberofGenes, rnd.Next(1, numberofGenes));
            int swappingLength = cutPoint_1.Length;
            int[] cutPoint_2 = new int[swappingLength];
            Array.Sort(cutPoint_1);
            int[] swappingPartFromF = new int[swappingLength];
            int[] swappingPartFromM = new int[swappingLength];
            int[] remainpartforA = new int[numberofGenes - swappingLength];
            int[] remainpartforB = new int[numberofGenes - swappingLength];

            int s;
            int r;
            int rf;
            int rm;

            s = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (cutPoint_1.Contains(i))
                {
                    swappingPartFromF[s] = Chromosomes[father][i];
                    s++;
                }
            }
            r = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (swappingPartFromF.Contains(Chromosomes[mother][i]) == true)
                {
                    cutPoint_2[r] = i;
                    r++;
                }

            }
            s = 0;
            for (int i = 0; i < cutPoint_2.Length; i++)
            {
                swappingPartFromM[s] = Chromosomes[mother][cutPoint_2[i]];
                s++;
            }


            rf = 0;
            rm = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (swappingPartFromF.Contains(Chromosomes[mother][i]) == false)
                {
                    remainpartforA[rf] = Chromosomes[mother][i];
                    rf++;
                }
                if (swappingPartFromM.Contains(Chromosomes[father][i]) == false)
                {
                    remainpartforB[rm] = Chromosomes[father][i];
                    rm++;
                }
            }

            rf = 0;
            rm = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (cutPoint_1.Contains(i) == true)
                {
                    Chromosomes[childb][i] = swappingPartFromM[rm];
                    rm++;
                }
                if (cutPoint_2.Contains(i) == true)
                {
                    Chromosomes[childa][i] = swappingPartFromF[rf];
                    rf++;
                }
            }

            rf = 0;
            rm = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (cutPoint_1.Contains(i) == false)
                {
                    Chromosomes[childb][i] = remainpartforB[rm];
                    rm++;
                }

                if (cutPoint_2.Contains(i) == false)
                {
                    Chromosomes[childa][i] = remainpartforA[rf];
                    rf++;
                }

            }


        }
        void CycleCrossover(int father, int mother, int childa, int childb)
        {
            int cutPoint = rnd.Next(1, numberofGenes);
            List<int> cutcycleelementofF = new List<int>();
            List<int> cutcycleindexofF = new List<int>();

            int elementofF = Chromosomes[father][cutPoint];
            int elementindexofF = Array.IndexOf(Chromosomes[father], elementofF);
            cutcycleelementofF.Add(elementofF);
            cutcycleindexofF.Add(elementindexofF);
            while (elementofF != -1)
            {
                elementindexofF = Array.IndexOf(Chromosomes[mother], cutcycleindexofF.Last());
                elementofF = Chromosomes[father][elementindexofF];
                if (cutcycleelementofF.Contains(elementofF))
                {
                    elementofF = -1;
                    break;
                }

                cutcycleelementofF.Add(elementofF);
                cutcycleindexofF.Add(elementindexofF);
            }
            int[] cycleSwappingIndex = cutcycleindexofF.ToArray();
            int swappingLength = cycleSwappingIndex.Length;
            int[] swappingPartfromF = new int[swappingLength];
            int[] swappingPartfromM = new int[swappingLength];
            int sindex = 0;
            foreach (var item in cycleSwappingIndex)
            {
                swappingPartfromF[sindex] = Chromosomes[father][item];
                swappingPartfromM[sindex] = Chromosomes[mother][item];
                sindex++;
            }

            int[] remainpartforA = new int[numberofGenes - swappingLength];
            int[] remainpartforB = new int[numberofGenes - swappingLength];
            int rindex = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (cycleSwappingIndex.Contains(i) == false)
                {
                    remainpartforA[rindex] = Chromosomes[mother][i];
                    remainpartforB[rindex] = Chromosomes[father][i];
                    rindex++;
                }
            }

            rindex = 0;
            sindex = 0;

            for (int i = 0; i < numberofGenes; i++)
            {
                if (cycleSwappingIndex.Contains(i) == true)
                {
                    Chromosomes[childa][i] = swappingPartfromF[sindex];
                    Chromosomes[childb][i] = swappingPartfromM[sindex];
                    sindex++;
                }
                else
                {
                    Chromosomes[childa][i] = remainpartforA[rindex];
                    Chromosomes[childb][i] = remainpartforB[rindex];
                    rindex++;
                }
            }
        }
        void InversionMutation(int beforemutation, int aftermutation, bool[]mutatedFlag)
        {
            int[] twoPoints = ShuffleaArray(numberofGenes, 2);
            Array.Sort(twoPoints);
            int p1 = twoPoints[0];
            int p2 = twoPoints[1];

            int inversedLength = (p2 - p1);
            int[] inversedPart = new int[inversedLength];
            for (int i = 0; i < inversedLength; i++)
            {
                inversedPart[i] = Chromosomes[beforemutation][p2 - i - 1];
            }
            int rindex = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if (i >= p1 && i < p2)
                {
                    Chromosomes[aftermutation][i] = inversedPart[rindex];
                    rindex++;
                }
                else
                {
                    Chromosomes[aftermutation][i] = Chromosomes[beforemutation][i];
                }
            }
        }

        void InsertionMuation(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            int mutationPoint = rnd.Next(0, numberofGenes);

            int[] remainingPart = new int[numberofGenes - 1];
            int rindex = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if(mutationPoint != Chromosomes[beforemutation][i])
                {
                    remainingPart[rindex] = Chromosomes[beforemutation][i];
                    rindex++;
                }
            }

            rindex = 0;
            int insertionPoint = rnd.Next(0, remainingPart.Length);
            for (int i = 0; i < numberofGenes; i++)
            {
                if(i != insertionPoint)
                {
                    Chromosomes[aftermutation][i] = remainingPart[rindex];
                    rindex++;
                }
                else
                {
                    Chromosomes[aftermutation][i] = mutationPoint;
                }
            }
        }

        void DisplacementMutation(int beforemutation, int aftermutation, bool[] mutationFlag)
        {
            int[] twoPoints = ShuffleaArray(numberofGenes, 2);
            Array.Sort(twoPoints);
            int p1 = twoPoints[0];
            int p2 = twoPoints[1];
            int displacementLength = (p2 - p1);

            int[] displacementPart = new int[displacementLength];
            for (int i = 0; i < displacementLength; i++)
            {
                displacementPart[i] = Chromosomes[beforemutation][p1 + i];
            }

            int[] remainingPart = new int[numberofGenes - displacementLength];
            int rindex = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if(displacementPart.Contains(Chromosomes[beforemutation][i]) == false)
                {
                    remainingPart[rindex] = Chromosomes[beforemutation][i];
                    rindex++;
                }
            }
            int insertPoint = rnd.Next(0, (numberofGenes - displacementLength));

            int dindex = 0;
            rindex = 0;
            for (int i = 0; i < numberofGenes; i++)
            {
                if(i == insertPoint)
                {
                    for (int j = i; j < i +  displacementLength; j++)
                    {
                        Chromosomes[aftermutation][j] = displacementPart[dindex];
                        dindex++;
                    }
                    i = i + displacementLength - 1;
                }
                else
                {
                    Chromosomes[aftermutation][i] = remainingPart[rindex];
                    rindex++;
                }

            }
        }
        void ReciprocalExchangeMutation(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            int[] twopoints = ShuffleaArray(numberofGenes, 2);
            Array.Sort(twopoints);
            int p1 = twopoints[0];
            int p2 = twopoints[1];

            for (int i = 0; i < numberofGenes; i++)
            {
                Chromosomes[aftermutation][i] = Chromosomes[beforemutation][i];
            }
            (Chromosomes[aftermutation][p2], Chromosomes[aftermutation][p1])
                = (Chromosomes[beforemutation][p1], Chromosomes[beforemutation][p2]);

        }

        void HeuristicMutation(int beforemutation, int aftermutation, bool[]mutatedFlag)
        {

        }
        public override void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb)
        {
            switch (CrossoverType)
            {
                case PermutationCrossoverType.ParialMappedX:
                    ParialMappedCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.OrderX:
                    OrderCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.PositionBasedX:
                    PostitionBasedCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.OrderBasedX:
                    OrderBasedCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.CycleX:
                    CycleCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.SubtourExchange:
                    break;
                default:
                    break;       
            }
            ObjectiveValue[childa] = objectivefunction(Chromosomes[childa]);
            ObjectiveValue[childb] = objectivefunction(Chromosomes[childb]);
        }


        public override void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            switch (PermuatationMutationType)
            {
                case PermuatationMutationType.Inversion:
                    InversionMutation(beforemutation, aftermutation, mutatedFlag);
                    break;
                case PermuatationMutationType.Insertion:
                    InsertionMuation(beforemutation, aftermutation, mutatedFlag);
                    break;
                case PermuatationMutationType.Displacement:
                    DisplacementMutation(beforemutation, aftermutation, mutatedFlag);
                    break;
                case PermuatationMutationType.ReciprocalExchange:
                    ReciprocalExchangeMutation(beforemutation, aftermutation, mutatedFlag);
                    break;
                case PermuatationMutationType.Heuristic:
                    HeuristicMutation(beforemutation, aftermutation, mutatedFlag);
                    break;
                default:
                    break;
            
            }
            ObjectiveValue[aftermutation] = objectivefunction(Chromosomes[aftermutation]);

        }


    }
}
