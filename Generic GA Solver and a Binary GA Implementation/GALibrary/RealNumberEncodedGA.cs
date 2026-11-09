using System;

namespace GALibrary
{
    /// <summary>Crossover operators for real-valued chromosomes.</summary>
    public enum RealNumberCrossoverType
    {
        /// <summary>Linear variation directed to the larger parent.</summary>
        LVD,

        /// <summary>Linear variation directed to the smaller parent.</summary>
        SVD,

        /// <summary>Multi-parent oriented extended segment crossover.</summary>
        MOES,

        /// <summary>Two-parent extended segment crossover.</summary>
        TES,

        /// <summary>Flat block / blended mean crossover.</summary>
        FBMS
    }

    /// <summary>Mutation operators for real-valued chromosomes.</summary>
    public enum RealNUmberMutationType
    {
        /// <summary>Non-uniform (dynamic) mutation whose step shrinks each generation.</summary>
        DynamicMutation
    }

    /// <summary>
    /// Real-valued encoding. Each gene is a double constrained to
    /// [lowerBound[i], upperBound[i]]. The class used to be internal and never
    /// overrode the crossover/mutation contract; it is now public and fully
    /// implements it.
    /// </summary>
    public class RealNumberEncodedGA : GenericSolver<double>
    {
        private readonly double[] lowerBound;
        private readonly double[] upperBound;
        private double degreeOfNonuniformity = 1.5;

        public RealNumberCrossoverType CrossoverType { get; set; } = RealNumberCrossoverType.LVD;

        /// <summary>Mutation operator used by this encoding.</summary>
        public RealNUmberMutationType RealMutationType { get; set; } = RealNUmberMutationType.DynamicMutation;

        /// <summary>Shape parameter of the non-uniform mutation (larger = stronger shrink).</summary>
        public double DegreeOfNonuniformity
        {
            get { return degreeOfNonuniformity; }
            set
            {
                if (double.IsNaN(value) || double.IsInfinity(value) || value <= 0.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "DegreeOfNonuniformity must be positive and finite.");
                }
                degreeOfNonuniformity = value;
            }
        }

        public double[] LowerBounds
        {
            get { return (double[])lowerBound.Clone(); }
        }

        public double[] UpperBounds
        {
            get { return (double[])upperBound.Clone(); }
        }

        public RealNumberEncodedGA(
            int numberofGenes,
            double[] lowerBound,
            double[] upperBound,
            GAOptimizationType optimizationType,
            ObjectiveFunction<double> objectiveFunction)
            : base(numberofGenes, optimizationType, objectiveFunction)
        {
            if (lowerBound == null) throw new ArgumentNullException(nameof(lowerBound));
            if (upperBound == null) throw new ArgumentNullException(nameof(upperBound));
            if (lowerBound.Length != numberofGenes || upperBound.Length != numberofGenes)
            {
                throw new ArgumentException("Both bound arrays must have one entry per gene.");
            }
            for (int i = 0; i < numberofGenes; i++)
            {
                if (double.IsNaN(lowerBound[i]) || double.IsNaN(upperBound[i])
                    || double.IsInfinity(lowerBound[i]) || double.IsInfinity(upperBound[i])
                    || lowerBound[i] >= upperBound[i])
                {
                    throw new ArgumentException("Each lower bound must be finite and strictly smaller than its upper bound.");
                }
            }
            this.lowerBound = (double[])lowerBound.Clone();
            this.upperBound = (double[])upperBound.Clone();
        }

        public override void InitializePopulation()
        {
            for (int row = 0; row < PopulationSize; row++)
            {
                for (int gene = 0; gene < numberofGenes; gene++)
                {
                    Chromosomes[row][gene] = lowerBound[gene] + rnd.NextDouble() * (upperBound[gene] - lowerBound[gene]);
                }
                ObjectiveValue[row] = Evaluate(Chromosomes[row]);
            }
        }

        private double Clamp(int gene, double value)
        {
            if (double.IsNaN(value)) return lowerBound[gene];
            if (value < lowerBound[gene]) return lowerBound[gene];
            if (value > upperBound[gene]) return upperBound[gene];
            return value;
        }

        // -----------------------------------------------------------------
        // Crossover operators (bound-aware variants of the original formulas)
        // -----------------------------------------------------------------

        public override void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb)
        {
            ValidateIndices(father, mother, childa, childb);
            double[] dad = Chromosomes[father];
            double[] mom = Chromosomes[mother];
            double[] childOne = Chromosomes[childa];
            double[] childTwo = Chromosomes[childb];

            for (int gene = 0; gene < numberofGenes; gene++)
            {
                double alpha = rnd.NextDouble();
                double small = Math.Min(dad[gene], mom[gene]);
                double large = Math.Max(dad[gene], mom[gene]);
                double low = lowerBound[gene];
                double high = upperBound[gene];

                switch (CrossoverType)
                {
                    case RealNumberCrossoverType.LVD:
                        childOne[gene] = Clamp(gene, alpha * low + (1.0 - alpha) * large);
                        childTwo[gene] = Clamp(gene, alpha * large + (1.0 - alpha) * high);
                        break;
                    case RealNumberCrossoverType.SVD:
                        childOne[gene] = Clamp(gene, alpha * low + (1.0 - alpha) * small);
                        childTwo[gene] = Clamp(gene, alpha * small + (1.0 - alpha) * high);
                        break;
                    case RealNumberCrossoverType.MOES:
                        {
                            double mix = rnd.NextDouble();
                            childOne[gene] = Clamp(gene, alpha * small + (1.0 - alpha) * large);
                            childTwo[gene] = mix > 0.5
                                ? Clamp(gene, alpha * large + (1.0 - alpha) * high)
                                : Clamp(gene, alpha * low + (1.0 - alpha) * small);
                        }
                        break;
                    case RealNumberCrossoverType.TES:
                        childOne[gene] = Clamp(gene, alpha * low + (1.0 - alpha) * small);
                        childTwo[gene] = Clamp(gene, alpha * large + (1.0 - alpha) * high);
                        break;
                    case RealNumberCrossoverType.FBMS:
                        childOne[gene] = Clamp(gene, alpha * small + (1.0 - alpha) * large);
                        childTwo[gene] = Clamp(gene, alpha * large + (1.0 - alpha) * small);
                        break;
                    default:
                        throw new ArgumentOutOfRangeException("CrossoverType", CrossoverType, "Unknown crossover operator.");
                }
            }

            ObjectiveValue[childa] = Evaluate(childOne);
            ObjectiveValue[childb] = Evaluate(childTwo);
        }

        // -----------------------------------------------------------------
        // Mutation
        // -----------------------------------------------------------------

        public override void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            if (mutatedFlag == null) throw new ArgumentNullException(nameof(mutatedFlag));
            if (mutatedFlag.Length != numberofGenes)
            {
                throw new ArgumentException("The mutation mask must match the chromosome length.", nameof(mutatedFlag));
            }
            ValidateIndices(beforemutation, aftermutation);

            double progress = IterationLimit > 0 ? Math.Min(1.0, (double)IterationID / IterationLimit) : 1.0;
            double shaping = Math.Pow(1.0 - progress, degreeOfNonuniformity);

            for (int gene = 0; gene < numberofGenes; gene++)
            {
                double value = Chromosomes[beforemutation][gene];
                if (mutatedFlag[gene])
                {
                    double random = rnd.NextDouble();
                    double delta = (1.0 - Math.Pow(random, shaping)) * (upperBound[gene] - lowerBound[gene]);
                    if (rnd.NextDouble() < 0.5)
                    {
                        value += delta;
                    }
                    else
                    {
                        value -= delta;
                    }
                    value = Clamp(gene, value);
                }
                Chromosomes[aftermutation][gene] = value;
            }

            ObjectiveValue[aftermutation] = Evaluate(Chromosomes[aftermutation]);
        }

        private void ValidateIndices(params int[] indicesToCheck)
        {
            for (int i = 0; i < indicesToCheck.Length; i++)
            {
                int index = indicesToCheck[i];
                if (index < 0 || index >= Chromosomes.Length)
                {
                    throw new ArgumentOutOfRangeException(nameof(indicesToCheck), "Chromosome index is outside the 3N buffer.");
                }
            }
        }
    }
}
