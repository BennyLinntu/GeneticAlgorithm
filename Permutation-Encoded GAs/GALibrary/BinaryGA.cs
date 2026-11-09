using System;
using System.Collections.Generic;

namespace GALibrary
{
    /// <summary>Cut-point crossover variants for binary chromosomes.</summary>
    public enum BinaryCrossoverType
    {
        OnePointCut,
        TwoPointCut,
        NPointCut
    }

    /// <summary>
    /// Binary encoding. A chromosome is a fixed-length string of 0/1 genes. The
    /// selection use case interprets the string as a row-major n x n matrix and
    /// penalises rows or columns that do not contain exactly one selected cell.
    /// </summary>
    public class BinaryGA : GenericSolver<byte>
    {
        private BinaryCrossoverType crossoverType = BinaryCrossoverType.OnePointCut;

        public BinaryCrossoverType CrossoverType
        {
            get { return crossoverType; }
            set
            {
                if (value != BinaryCrossoverType.OnePointCut
                    && value != BinaryCrossoverType.TwoPointCut
                    && value != BinaryCrossoverType.NPointCut)
                {
                    throw new ArgumentOutOfRangeException(nameof(value));
                }
                crossoverType = value;
            }
        }

        public BinaryGA(
            int numberofGenes,
            GAOptimizationType optimizationType,
            ObjectiveFunction<byte> objectiveFunction,
            BinaryCrossoverType crossoverType)
            : base(numberofGenes, optimizationType, objectiveFunction)
        {
            CrossoverType = crossoverType;
        }

        // -----------------------------------------------------------------
        // Initialisation
        // -----------------------------------------------------------------

        public override void InitializePopulation()
        {
            for (int row = 0; row < PopulationSize; row++)
            {
                for (int gene = 0; gene < numberofGenes; gene++)
                {
                    Chromosomes[row][gene] = (byte)rnd.Next(2);
                }
                ObjectiveValue[row] = Evaluate(Chromosomes[row]);
            }
        }

        // -----------------------------------------------------------------
        // Constraints
        // -----------------------------------------------------------------

        /// <summary>
        /// Returns row violations followed by column violations for the square
        /// cost matrix represented by <paramref name="chrom"/>. If the
        /// chromosome length is not a perfect square there is no implied matrix
        /// and an empty array is returned.
        /// </summary>
        public override int[] ReturnChromosomesViolation(byte[] chrom)
        {
            if (chrom == null || chrom.Length != numberofGenes)
            {
                return new int[0];
            }

            int size = (int)Math.Round(Math.Sqrt(numberofGenes));
            if (size * size != numberofGenes)
            {
                return new int[0];
            }

            int[] violations = new int[2 * size];
            for (int row = 0; row < size; row++)
            {
                int count = 0;
                for (int column = 0; column < size; column++)
                {
                    count += chrom[row * size + column] == 0 ? 0 : 1;
                }
                violations[row] = Math.Abs(count - 1);
            }

            for (int column = 0; column < size; column++)
            {
                int count = 0;
                for (int row = 0; row < size; row++)
                {
                    count += chrom[row * size + column] == 0 ? 0 : 1;
                }
                violations[size + column] = Math.Abs(count - 1);
            }

            return violations;
        }

        // -----------------------------------------------------------------
        // Mutation
        // -----------------------------------------------------------------

        public override void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            if (mutatedFlag == null)
            {
                throw new ArgumentNullException(nameof(mutatedFlag));
            }
            if (mutatedFlag.Length != numberofGenes)
            {
                throw new ArgumentException("The mutation mask must match the chromosome length.", nameof(mutatedFlag));
            }
            if (beforemutation < 0 || beforemutation >= Chromosomes.Length
                || aftermutation < 0 || aftermutation >= Chromosomes.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(aftermutation), "Chromosome index is outside the 3N buffer.");
            }

            for (int gene = 0; gene < numberofGenes; gene++)
            {
                byte value = Chromosomes[beforemutation][gene];
                Chromosomes[aftermutation][gene] = value == 0 ? (byte)0 : (byte)1;
                if (mutatedFlag[gene])
                {
                    Chromosomes[aftermutation][gene] = (byte)(Chromosomes[aftermutation][gene] == 0 ? 1 : 0);
                }
            }
            ObjectiveValue[aftermutation] = Evaluate(Chromosomes[aftermutation]);
        }

        // -----------------------------------------------------------------
        // Crossover
        // -----------------------------------------------------------------

        public override void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb)
        {
            if (father < 0 || father >= Chromosomes.Length
                || mother < 0 || mother >= Chromosomes.Length
                || childa < 0 || childa >= Chromosomes.Length
                || childb < 0 || childb >= Chromosomes.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(childa), "Chromosome index is outside the 3N buffer.");
            }

            byte[] dad = Chromosomes[father];
            byte[] mom = Chromosomes[mother];
            byte[] childOne = Chromosomes[childa];
            byte[] childTwo = Chromosomes[childb];

            if (numberofGenes < 2)
            {
                Array.Copy(dad, childOne, numberofGenes);
                Array.Copy(mom, childTwo, numberofGenes);
                ObjectiveValue[childa] = Evaluate(childOne);
                ObjectiveValue[childb] = Evaluate(childTwo);
                return;
            }

            switch (crossoverType)
            {
                case BinaryCrossoverType.OnePointCut:
                    {
                        // A cut in [1, n-1] guarantees that both children receive
                        // genes from both parents (no degenerate empty segment).
                        int cut = rnd.Next(1, numberofGenes);
                        for (int gene = 0; gene < numberofGenes; gene++)
                        {
                            if (gene < cut)
                            {
                                childOne[gene] = dad[gene];
                                childTwo[gene] = mom[gene];
                            }
                            else
                            {
                                childOne[gene] = mom[gene];
                                childTwo[gene] = dad[gene];
                            }
                        }
                    }
                    break;

                case BinaryCrossoverType.TwoPointCut:
                    {
                        // Distinct, ordered cuts: 0 <= first < second <= n-1.
                        int first = rnd.Next(0, numberofGenes - 1);
                        int second = rnd.Next(first + 1, numberofGenes);
                        for (int gene = 0; gene < numberofGenes; gene++)
                        {
                            bool inside = gene >= first && gene < second;
                            if (inside)
                            {
                                childOne[gene] = dad[gene];
                                childTwo[gene] = mom[gene];
                            }
                            else
                            {
                                childOne[gene] = mom[gene];
                                childTwo[gene] = dad[gene];
                            }
                        }
                    }
                    break;

                case BinaryCrossoverType.NPointCut:
                    {
                        // At least one cut and at most n-1 cuts, all distinct and
                        // strictly inside the chromosome, so no child is a plain copy.
                        int cuts = rnd.Next(1, numberofGenes);
                        List<int> positions = new List<int>();
                        for (int i = 1; i < numberofGenes; i++)
                        {
                            positions.Add(i);
                        }
                        for (int i = positions.Count - 1; i > 0; i--)
                        {
                            int j = rnd.Next(i + 1);
                            int swap = positions[i];
                            positions[i] = positions[j];
                            positions[j] = swap;
                        }

                        int cutCount = Math.Min(cuts, positions.Count);
                        int[] cutPoints = new int[cutCount];
                        for (int i = 0; i < cutCount; i++)
                        {
                            cutPoints[i] = positions[i];
                        }
                        Array.Sort(cutPoints);

                        bool takeFather = true;
                        int cutIndex = 0;
                        for (int gene = 0; gene < numberofGenes; gene++)
                        {
                            while (cutIndex < cutPoints.Length && gene == cutPoints[cutIndex])
                            {
                                takeFather = !takeFather;
                                cutIndex++;
                            }
                            if (takeFather)
                            {
                                childOne[gene] = dad[gene];
                                childTwo[gene] = mom[gene];
                            }
                            else
                            {
                                childOne[gene] = mom[gene];
                                childTwo[gene] = dad[gene];
                            }
                        }
                    }
                    break;

                default:
                    throw new ArgumentOutOfRangeException("CrossoverType", crossoverType, "Unknown crossover operator.");
            }

            ObjectiveValue[childa] = Evaluate(childOne);
            ObjectiveValue[childb] = Evaluate(childTwo);
        }
    }
}

