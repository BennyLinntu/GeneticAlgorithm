using System;
using System.Drawing;
using System.Globalization;
using System.Text;
using System.Windows.Forms.DataVisualization.Charting;

namespace GALibrary
{
    // ---------------------------------------------------------------------
    // Shared enums and delegates used by every solver in the library.
    // ---------------------------------------------------------------------

    /// <summary>Direction of the optimisation.</summary>
    public enum GAOptimizationType
    {
        Minimization,
        Maximization
    }

    /// <summary>How mutation is counted.</summary>
    public enum GAMutationType
    {
        /// <summary>Mutation probability is applied per gene.</summary>
        GeneNumberBased,

        /// <summary>Mutation probability is applied per chromosome.</summary>
        ChromosomesNumberBased
    }

    /// <summary>Survivor / parent selection strategy.</summary>
    public enum GASelectionType
    {
        /// <summary>Truncation: keep the best <see cref="GenericSolver{T}.PopulationSize"/> individuals.</summary>
        Deterministic,

        /// <summary>Fitness-proportionate roulette-wheel selection.</summary>
        Stochastic
    }

    /// <summary>Legacy delegate kept for binary compatibility with the original selection code.</summary>
    public delegate double COPObjectiveFunction(double[] solution);

    /// <summary>
    /// Template-method genetic algorithm engine shared by the binary, permutation
    /// and real-valued encodings.
    ///
    /// Buffer layout (capacity = 3 * <see cref="PopulationSize"/>):
    /// <code>
    ///   [0 .. N)              current parents
    ///   [N .. N + C)          crossover children (C = NumberofCrossedChildren)
    ///   [N + C .. N + C + M)  mutation children  (M = NumberofMutatedChildren)
    /// </code>
    /// Every row is fully allocated in <see cref="Reset"/> and every operator is
    /// required to overwrite the complete child row, so no uninitialised genes
    /// can leak into an evaluation.
    /// </summary>
    public abstract class GenericSolver<T>
    {
        // -----------------------------------------------------------------
        // Random source. Public to preserve the original API; use
        // SetRandomSeed for reproducible runs.
        // -----------------------------------------------------------------
        public Random rnd = new Random();

        protected T[][] chromosomes;
        protected double[] objectiveValue;
        protected double[] fitnessValue;

        private T[] soFartheBestSolution;
        private double soFartheBestObjectValue;
        private double initialBestObjective;
        private double iterationBestObjective;
        private double iterationAverageObjective;

        protected int[] indices;
        protected T[][] selectedChromosome;
        protected double[] selectedObjectiveValue;
        protected int numberofGenes;

        protected ObjectiveFunction<T> objectivefunction;
        protected double lessFitnessFraction;

        // Parameters with validation in their setters.
        private int populationSize = 10;
        private double crossoverRate = 0.8;
        private double mutationRate = 0.2;
        private double leastFitnessFraction = 0.1;
        private double penaltyFactor = 100;
        private int iterationLimit = 200;

        private GAOptimizationType optimizationType = GAOptimizationType.Minimization;
        private GAMutationType mutationType = GAMutationType.ChromosomesNumberBased;
        private GASelectionType selectionType = GASelectionType.Deterministic;

        /// <summary>Objective function supplied by the client.</summary>
        public delegate double ObjectiveFunction<S>(S[] solution);

        private int numberofCrossedChildren;
        private int numberofMutatedChildren;

        private Series seriesSoFartheBest;
        private Series seriesAverage;
        private Series seriesIterationtheBest;

        private bool[][] mutatedGenes;
        private int iterationID;

        // -----------------------------------------------------------------
        // Construction
        // -----------------------------------------------------------------

        protected GenericSolver(int numberofGenes, GAOptimizationType optimizationType, ObjectiveFunction<T> objectiveFunction)
        {
            if (numberofGenes < 1)
            {
                throw new ArgumentOutOfRangeException(nameof(numberofGenes), "A chromosome must contain at least one gene.");
            }
            if (objectiveFunction == null)
            {
                throw new ArgumentNullException(nameof(objectiveFunction));
            }
            if (optimizationType != GAOptimizationType.Minimization && optimizationType != GAOptimizationType.Maximization)
            {
                throw new ArgumentOutOfRangeException(nameof(optimizationType));
            }

            this.numberofGenes = numberofGenes;
            this.optimizationType = optimizationType;
            this.objectivefunction = objectiveFunction;

            seriesSoFartheBest = new Series("So Far the Best");
            seriesSoFartheBest.ChartType = SeriesChartType.Line;
            seriesSoFartheBest.Color = Color.Red;
            seriesSoFartheBest.BorderWidth = 2;
            seriesSoFartheBest.IsVisibleInLegend = true;

            seriesAverage = new Series("Population Average");
            seriesAverage.ChartType = SeriesChartType.Line;
            seriesAverage.Color = Color.Green;
            seriesAverage.BorderWidth = 2;
            seriesAverage.IsVisibleInLegend = true;

            seriesIterationtheBest = new Series("Iteration Best");
            seriesIterationtheBest.ChartType = SeriesChartType.Line;
            seriesIterationtheBest.Color = Color.Blue;
            seriesIterationtheBest.BorderWidth = 2;
            seriesIterationtheBest.IsVisibleInLegend = true;
        }

        // -----------------------------------------------------------------
        // Validated public properties (names preserved from the original API)
        // -----------------------------------------------------------------

        public int PopulationSize
        {
            get { return populationSize; }
            set
            {
                if (value < 2)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "PopulationSize must be at least 2.");
                }
                if (value > 100000)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "PopulationSize must not exceed 100000.");
                }
                populationSize = value;
            }
        }

        public double CrossoverRate
        {
            get { return crossoverRate; }
            set
            {
                if (double.IsNaN(value) || value < 0.0 || value > 1.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "CrossoverRate must be in [0, 1].");
                }
                crossoverRate = value;
            }
        }

        public double MutationRate
        {
            get { return mutationRate; }
            set
            {
                if (double.IsNaN(value) || value < 0.0 || value > 1.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "MutationRate must be in [0, 1].");
                }
                mutationRate = value;
            }
        }

        public T[][] Chromosomes
        {
            get { return chromosomes; }
            set { chromosomes = value; }
        }

        public GASelectionType SelectionType
        {
            get { return selectionType; }
            set { selectionType = value; }
        }

        /// <summary>Best solution seen since the last <see cref="Reset"/>.</summary>
        public T[] SoFartheBestSoltion
        {
            get { return soFartheBestSolution; }
            set { soFartheBestSolution = value; }
        }

        /// <summary>Correctly spelled alias of <see cref="SoFartheBestSoltion"/>.</summary>
        public T[] SoFarBestSolution
        {
            get { return soFartheBestSolution; }
            set { soFartheBestSolution = value; }
        }

        public double SoFartheBestObjectSolutionValue
        {
            get { return soFartheBestObjectValue; }
            set { soFartheBestObjectValue = value; }
        }

        public int IterationID
        {
            get { return iterationID; }
            set { iterationID = value; }
        }

        public int NumberofCrossedChildren
        {
            get { return numberofCrossedChildren; }
            set { numberofCrossedChildren = value; }
        }

        public int NumberofMutatedChildren
        {
            get { return numberofMutatedChildren; }
            set { numberofMutatedChildren = value; }
        }

        public Series SeriesSoFartheBest
        {
            get { return seriesSoFartheBest; }
            set { seriesSoFartheBest = value; }
        }

        public Series SeriesAverage
        {
            get { return seriesAverage; }
            set { seriesAverage = value; }
        }

        public Series SeriesIterationtheBest
        {
            get { return seriesIterationtheBest; }
            set { seriesIterationtheBest = value; }
        }

        public int IterationLimit
        {
            get { return iterationLimit; }
            set
            {
                if (value < 0)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "IterationLimit must not be negative.");
                }
                if (value > 1000000)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "IterationLimit must not exceed 1000000.");
                }
                iterationLimit = value;
            }
        }

        public double PenaltyFactor
        {
            get { return penaltyFactor; }
            set
            {
                if (double.IsNaN(value) || double.IsInfinity(value) || value < 0.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "PenaltyFactor must be a finite, non-negative number.");
                }
                penaltyFactor = value;
            }
        }

        /// <summary>Legacy (misspelled) property retained for the original UI.</summary>
        public double LeastFinessFraction
        {
            get { return leastFitnessFraction; }
            set
            {
                if (double.IsNaN(value) || value <= 0.0 || value > 1.0)
                {
                    throw new ArgumentOutOfRangeException(nameof(value), "Least fitness fraction must be in (0, 1].");
                }
                leastFitnessFraction = value;
            }
        }

        /// <summary>Correctly spelled alias of <see cref="LeastFinessFraction"/>.</summary>
        public double LeastFitnessFraction
        {
            get { return leastFitnessFraction; }
            set { LeastFinessFraction = value; }
        }

        public GAOptimizationType OptimizationType
        {
            get { return optimizationType; }
            set
            {
                if (value != GAOptimizationType.Minimization && value != GAOptimizationType.Maximization)
                {
                    throw new ArgumentOutOfRangeException(nameof(value));
                }
                optimizationType = value;
            }
        }

        public GAMutationType MutationType
        {
            get { return mutationType; }
            set { mutationType = value; }
        }

        /// <summary>Legacy alias of <see cref="ObjectiveValue"/>.</summary>
        public double[] ObjectiveValue1
        {
            get { return objectiveValue; }
            set { objectiveValue = value; }
        }

        public double[] ObjectiveValue
        {
            get { return objectiveValue; }
            set { objectiveValue = value; }
        }

        public double InitialBestObjective
        {
            get { return initialBestObjective; }
            set { initialBestObjective = value; }
        }

        /// <summary>Best objective in the current parent + children pool.</summary>
        public double IterationBestObjective
        {
            get { return iterationBestObjective; }
            set { iterationBestObjective = value; }
        }

        /// <summary>Average objective in the current parent + children pool.</summary>
        public double IterationAverageObjective
        {
            get { return iterationAverageObjective; }
            set { iterationAverageObjective = value; }
        }

        /// <summary>Number of rows currently in use (N + C + M).</summary>
        public int ActivePoolSize
        {
            get { return populationSize + numberofCrossedChildren + numberofMutatedChildren; }
        }

        public int NumberofGenes
        {
            get { return numberofGenes; }
        }

        // -----------------------------------------------------------------
        // Objective evaluation, fitness scaling and constraint penalties
        // -----------------------------------------------------------------

        /// <summary>Raw objective value with no constraint penalty.</summary>
        public double ComputeRawObjective(T[] chromosome)
        {
            if (chromosome == null)
            {
                throw new ArgumentNullException(nameof(chromosome));
            }
            double value = objectivefunction(chromosome);
            if (double.IsNaN(value))
            {
                return optimizationType == GAOptimizationType.Maximization ? double.MinValue / 2.0 : double.MaxValue / 2.0;
            }
            if (double.IsPositiveInfinity(value))
            {
                return optimizationType == GAOptimizationType.Maximization ? double.MaxValue / 2.0 : double.MaxValue / 2.0;
            }
            if (double.IsNegativeInfinity(value))
            {
                return optimizationType == GAOptimizationType.Maximization ? double.MinValue / 2.0 : double.MinValue / 2.0;
            }
            return value;
        }

        /// <summary>
        /// Objective value plus (minimisation) or minus (maximisation) the weighted
        /// constraint violation reported by <see cref="ReturnChromosomesViolation"/>.
        /// </summary>
        public double ComputePenalizedObjective(T[] chromosome)
        {
            double value = ComputeRawObjective(chromosome);
            if (penaltyFactor > 0.0)
            {
                int[] violations = ReturnChromosomesViolation(chromosome);
                if (violations != null && violations.Length > 0)
                {
                    double total = 0.0;
                    for (int i = 0; i < violations.Length; i++)
                    {
                        total += Math.Abs(violations[i]);
                    }
                    if (total > 0.0)
                    {
                        value += optimizationType == GAOptimizationType.Minimization
                            ? penaltyFactor * total
                            : -penaltyFactor * total;
                    }
                }
            }
            return value;
        }

        /// <summary>Evaluates a chromosome and applies the constraint penalty.</summary>
        protected double Evaluate(T[] chromosome)
        {
            return ComputePenalizedObjective(chromosome);
        }

        /// <summary>
        /// Converts objective values into strictly positive fitness values.
        /// Minimisation ranks the smallest objective highest; maximisation ranks
        /// the largest objective highest. The "least fitness fraction" ensures
        /// that the worst individual still keeps a small chance of being sampled.
        /// </summary>
        public void SetFitnessandObjectives(int total, double leastFitnessFraction)
        {
            EnsureReady();
            if (total < 1 || total > chromosomes.Length)
            {
                throw new ArgumentOutOfRangeException(nameof(total));
            }
            if (double.IsNaN(leastFitnessFraction) || leastFitnessFraction <= 0.0 || leastFitnessFraction > 1.0)
            {
                throw new ArgumentOutOfRangeException(nameof(leastFitnessFraction));
            }
            if (fitnessValue == null || fitnessValue.Length < total)
            {
                fitnessValue = new double[chromosomes.Length];
            }

            double objectiveMin = double.MaxValue;
            double objectiveMax = double.MinValue;
            for (int i = 0; i < total; i++)
            {
                double value = objectiveValue[i];
                if (double.IsNaN(value) || double.IsInfinity(value))
                {
                    value = optimizationType == GAOptimizationType.Maximization ? double.MinValue / 2.0 : double.MaxValue / 2.0;
                    objectiveValue[i] = value;
                }
                if (value < objectiveMin) objectiveMin = value;
                if (value > objectiveMax) objectiveMax = value;
            }

            double range = objectiveMax - objectiveMin;
            if (double.IsNaN(range) || double.IsInfinity(range))
            {
                range = 0.0;
            }
            double beta = Math.Max(leastFitnessFraction * range, 1e-9);

            if (optimizationType == GAOptimizationType.Maximization)
            {
                for (int i = 0; i < total; i++)
                {
                    fitnessValue[i] = beta + (objectiveValue[i] - objectiveMin);
                }
            }
            else
            {
                for (int i = 0; i < total; i++)
                {
                    fitnessValue[i] = beta + (objectiveMax - objectiveValue[i]);
                }
            }
        }

        // -----------------------------------------------------------------
        // Selection
        // -----------------------------------------------------------------

        /// <summary>
        /// Selects <see cref="PopulationSize"/> survivors from the combined
        /// parent + crossover + mutation pool.
        /// </summary>
        public void PerformSelectionOperation()
        {
            EnsureReady();
            int total = ActivePoolSize;
            if (total < 1)
            {
                throw new InvalidOperationException("The pool is empty.");
            }

            SetFitnessandObjectives(total, leastFitnessFraction);

            // Sort by fitness descending; ties keep the lower index first so the
            // result is deterministic for a fixed seed.
            int[] order = new int[total];
            for (int i = 0; i < total; i++)
            {
                order[i] = i;
            }
            Array.Sort(order, delegate (int a, int b)
            {
                int byFitness = fitnessValue[b].CompareTo(fitnessValue[a]);
                if (byFitness != 0) return byFitness;
                return a.CompareTo(b);
            });

            if (selectionType == GASelectionType.Stochastic)
            {
                StochasticSelection(order, total);
            }
            else
            {
                for (int i = 0; i < populationSize; i++)
                {
                    CopyToSelection(order[i], i);
                }
            }
        }

        private void StochasticSelection(int[] order, int total)
        {
            // Roulette wheel over the (already sorted) fitness values. The wheel
            // is cumulative and each of the N draws uses an independent dart, so
            // fitter individuals occupy larger slices of the wheel.
            double sum = 0.0;
            for (int i = 0; i < total; i++)
            {
                double f = fitnessValue[i];
                if (double.IsNaN(f) || double.IsInfinity(f) || f < 0.0)
                {
                    f = 0.0;
                }
                sum += f;
            }
            if (!(sum > 0.0) || double.IsInfinity(sum))
            {
                // Degenerate wheel: fall back to truncation.
                for (int i = 0; i < populationSize; i++)
                {
                    CopyToSelection(order[i], i);
                }
                return;
            }

            double[] cumulative = new double[total];
            double running = 0.0;
            for (int i = 0; i < total; i++)
            {
                double f = fitnessValue[i];
                if (double.IsNaN(f) || double.IsInfinity(f) || f < 0.0) f = 0.0;
                running += f;
                cumulative[i] = running;
            }

            for (int slot = 0; slot < populationSize; slot++)
            {
                double dart = rnd.NextDouble() * running;
                int pick = total - 1;
                for (int i = 0; i < total; i++)
                {
                    if (dart < cumulative[i])
                    {
                        pick = i;
                        break;
                    }
                }
                CopyToSelection(pick, slot);
            }
        }

        private void CopyToSelection(int source, int destination)
        {
            Array.Copy(chromosomes[source], selectedChromosome[destination], numberofGenes);
            selectedObjectiveValue[destination] = objectiveValue[source];
        }

        private void CopyAllSelectionChromosomes()
        {
            for (int i = 0; i < populationSize; i++)
            {
                Array.Copy(selectedChromosome[i], chromosomes[i], numberofGenes);
                objectiveValue[i] = selectedObjectiveValue[i];
            }
        }

        // -----------------------------------------------------------------
        // Best-so-far tracking
        // -----------------------------------------------------------------

        private void UpdateSoFartheBestandIterationtheBest()
        {
            int total = ActivePoolSize;
            int bestIndex = 0;
            for (int i = 1; i < total; i++)
            {
                if (IsBetter(objectiveValue[i], objectiveValue[bestIndex]))
                {
                    bestIndex = i;
                }
            }

            iterationBestObjective = objectiveValue[bestIndex];

            double average = 0.0;
            int count = 0;
            for (int i = 0; i < total; i++)
            {
                double value = objectiveValue[i];
                if (!double.IsNaN(value) && !double.IsInfinity(value))
                {
                    average += value;
                    count++;
                }
            }
            iterationAverageObjective = count > 0 ? average / count : double.NaN;

            bool unset = double.IsNaN(soFartheBestObjectValue)
                || (optimizationType == GAOptimizationType.Minimization && soFartheBestObjectValue == double.MaxValue)
                || (optimizationType == GAOptimizationType.Maximization && soFartheBestObjectValue == double.MinValue);

            if (unset || IsBetter(objectiveValue[bestIndex], soFartheBestObjectValue))
            {
                soFartheBestObjectValue = objectiveValue[bestIndex];
                if (soFartheBestSolution == null || soFartheBestSolution.Length != numberofGenes)
                {
                    soFartheBestSolution = new T[numberofGenes];
                }
                Array.Copy(chromosomes[bestIndex], soFartheBestSolution, numberofGenes);
            }
        }

        private bool IsBetter(double candidate, double reference)
        {
            if (double.IsNaN(candidate)) return false;
            if (double.IsNaN(reference)) return true;
            return optimizationType == GAOptimizationType.Maximization
                ? candidate > reference
                : candidate < reference;
        }

        // -----------------------------------------------------------------
        // Crossover / mutation drivers
        // -----------------------------------------------------------------

        protected void PerformCrossoverOperation()
        {
            EnsureReady();
            int requested = (int)(populationSize * crossoverRate);
            if (requested % 2 == 1)
            {
                requested--;
            }
            if (requested > populationSize)
            {
                requested = populationSize - (populationSize % 2);
            }
            if (requested < 0)
            {
                requested = 0;
            }
            NumberofCrossedChildren = requested;

            if (NumberofCrossedChildren < 2)
            {
                NumberofCrossedChildren = 0;
                return;
            }

            int pairs = NumberofCrossedChildren / 2;
            int[] shuffled = ShuffleIndiceArray(populationSize);
            for (int p = 0; p < pairs; p++)
            {
                int father = shuffled[2 * p];
                int mother = shuffled[2 * p + 1];
                int childA = populationSize + 2 * p;
                int childB = populationSize + 2 * p + 1;
                GenerateCrossoveredChromosomes(father, mother, childA, childB);
            }
        }

        protected void PerformMutationOperation()
        {
            EnsureReady();
            if (mutatedGenes == null || mutatedGenes.Length != populationSize)
            {
                mutatedGenes = new bool[populationSize][];
                for (int i = 0; i < populationSize; i++)
                {
                    mutatedGenes[i] = new bool[numberofGenes];
                }
            }
            for (int row = 0; row < populationSize; row++)
            {
                for (int gene = 0; gene < numberofGenes; gene++)
                {
                    mutatedGenes[row][gene] = false;
                }
            }

            // The 3N buffer holds N parents, C crossover children and M mutation
            // children, therefore M <= 2N - C. In particular, mutation can still
            // create one child per parent when crossover used half of the buffer.
            int capacity = 2 * populationSize - NumberofCrossedChildren;
            if (capacity < 0) capacity = 0;
            if (capacity > populationSize) capacity = populationSize;

            if (mutationType == GAMutationType.GeneNumberBased)
            {
                int genePool = populationSize * numberofGenes;
                int mutatedGeneCount = (int)Math.Round(mutationRate * genePool, MidpointRounding.AwayFromZero);
                if (mutationRate > 0.0 && mutatedGeneCount < 1)
                {
                    mutatedGeneCount = 1;
                }
                if (mutatedGeneCount > genePool) mutatedGeneCount = genePool;

                int[] slots = ShuffleIndiceArray(genePool);
                for (int k = 0; k < mutatedGeneCount; k++)
                {
                    int row = slots[k] / numberofGenes;
                    int gene = slots[k] % numberofGenes;
                    mutatedGenes[row][gene] = true;
                }

                int created = 0;
                for (int row = 0; row < populationSize && created < capacity; row++)
                {
                    bool any = false;
                    for (int gene = 0; gene < numberofGenes; gene++)
                    {
                        if (mutatedGenes[row][gene])
                        {
                            any = true;
                            break;
                        }
                    }
                    if (!any) continue;

                    int destination = populationSize + NumberofCrossedChildren + created;
                    GenerateMutatedChromosomes(row, destination, mutatedGenes[row]);
                    created++;
                }
                NumberofMutatedChildren = created;
            }
            else
            {
                int mutatedChromosomeCount = (int)Math.Round(mutationRate * populationSize, MidpointRounding.AwayFromZero);
                if (mutationRate > 0.0 && mutatedChromosomeCount < 1)
                {
                    mutatedChromosomeCount = 1;
                }
                if (mutatedChromosomeCount > capacity) mutatedChromosomeCount = capacity;

                int[] shuffled = ShuffleIndiceArray(populationSize);
                for (int i = 0; i < mutatedChromosomeCount; i++)
                {
                    int row = shuffled[i];
                    int gene = rnd.Next(numberofGenes);
                    mutatedGenes[row][gene] = true;
                    int destination = populationSize + NumberofCrossedChildren + i;
                    GenerateMutatedChromosomes(row, destination, mutatedGenes[row]);
                }
                NumberofMutatedChildren = mutatedChromosomeCount;
            }
        }

        // -----------------------------------------------------------------
        // Contract implemented by the encodings
        // -----------------------------------------------------------------

        /// <summary>Creates the initial population and evaluates it.</summary>
        public abstract void InitializePopulation();

        /// <summary>Fills both child rows with a recombination of two parents.</summary>
        public abstract void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb);

        /// <summary>Fills the child row with a mutated copy of the parent row.</summary>
        public abstract void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag);

        /// <summary>
        /// Constraint violations of a chromosome. Encodings without constraints
        /// return an empty array. The binary encoding returns row and column
        /// violations of the cost matrix.
        /// </summary>
        public virtual int[] ReturnChromosomesViolation(T[] chrom)
        {
            return new int[0];
        }

        /// <summary>Resets the internal buffers. Derived encodings may clear extra state.</summary>
        public virtual void ClearChromsomesandObjectiveArray()
        {
            if (chromosomes == null) return;
            for (int row = 0; row < chromosomes.Length; row++)
            {
                if (chromosomes[row] != null)
                {
                    Array.Clear(chromosomes[row], 0, chromosomes[row].Length);
                }
                if (objectiveValue != null && row < objectiveValue.Length)
                {
                    objectiveValue[row] = 0.0;
                }
            }
        }

        // -----------------------------------------------------------------
        // Lifecycle
        // -----------------------------------------------------------------

        /// <summary>Allocates the 3N buffers, seeds the population and prepares generation 0.</summary>
        public void Reset()
        {
            if (objectivefunction == null)
            {
                throw new InvalidOperationException("No objective function has been supplied.");
            }
            if (numberofGenes < 1)
            {
                throw new InvalidOperationException("The chromosome length must be positive.");
            }

            iterationID = 0;
            numberofCrossedChildren = 0;
            numberofMutatedChildren = 0;

            int threeTimesSize = populationSize * 3;
            indices = new int[threeTimesSize];
            objectiveValue = new double[threeTimesSize];
            fitnessValue = new double[threeTimesSize];
            chromosomes = new T[threeTimesSize][];
            for (int row = 0; row < threeTimesSize; row++)
            {
                chromosomes[row] = new T[numberofGenes];
            }

            selectedObjectiveValue = new double[populationSize];
            selectedChromosome = new T[populationSize][];
            for (int row = 0; row < populationSize; row++)
            {
                selectedChromosome[row] = new T[numberofGenes];
            }

            mutatedGenes = new bool[populationSize][];
            for (int row = 0; row < populationSize; row++)
            {
                mutatedGenes[row] = new bool[numberofGenes];
            }

            InitializePopulation();

            soFartheBestSolution = new T[numberofGenes];
            soFartheBestObjectValue = optimizationType == GAOptimizationType.Maximization
                ? double.MinValue
                : double.MaxValue;
            iterationBestObjective = double.NaN;
            iterationAverageObjective = double.NaN;

            UpdateSoFartheBestandIterationtheBest();
            initialBestObjective = soFartheBestObjectValue;

            PerformSelectionOperation();

            seriesAverage.Points.Clear();
            seriesIterationtheBest.Points.Clear();
            seriesSoFartheBest.Points.Clear();
            SeriesUpdate();
        }

        /// <summary>Runs a single generation, respecting <see cref="IterationLimit"/>.</summary>
        public void RunOneIteration()
        {
            EnsureReady();
            if (iterationID >= iterationLimit)
            {
                return;
            }

            if (iterationID != 0)
            {
                CopyAllSelectionChromosomes();
            }

            PerformCrossoverOperation();
            PerformMutationOperation();
            UpdateSoFartheBestandIterationtheBest();
            PerformSelectionOperation();

            iterationID++;
            SeriesUpdate();
        }

        /// <summary>
        /// Runs at most <paramref name="iteration"/> generations but never more
        /// than the remaining budget defined by <see cref="IterationLimit"/>.
        /// </summary>
        public void RunToEnd(int iteration)
        {
            EnsureReady();
            if (iteration < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(iteration));
            }
            long requestedEnd = (long)iterationID + iteration;
            long end = Math.Min((long)iterationLimit, requestedEnd);
            while (iterationID < end)
            {
                RunOneIteration();
            }
        }

        private void SeriesUpdate()
        {
            if (seriesAverage == null || seriesIterationtheBest == null || seriesSoFartheBest == null)
            {
                return;
            }
            if (seriesAverage.Points.Count > 100000)
            {
                seriesAverage.Points.Clear();
                seriesIterationtheBest.Points.Clear();
                seriesSoFartheBest.Points.Clear();
            }
            if (!double.IsNaN(iterationAverageObjective))
            {
                seriesAverage.Points.AddXY(iterationID, iterationAverageObjective);
            }
            if (!double.IsNaN(iterationBestObjective))
            {
                seriesIterationtheBest.Points.AddXY(iterationID, iterationBestObjective);
            }
            if (!double.IsNaN(soFartheBestObjectValue))
            {
                seriesSoFartheBest.Points.AddXY(iterationID, soFartheBestObjectValue);
            }
        }

        // -----------------------------------------------------------------
        // Helpers
        // -----------------------------------------------------------------

        /// <summary>Seeds the random source so that a run can be reproduced.</summary>
        public void SetRandomSeed(int seed)
        {
            rnd = new Random(seed);
        }

        /// <summary>Returns a Fisher-Yates shuffle of 0..limit-1 using the solver RNG.</summary>
        public int[] ShuffleIndiceArray(int limit)
        {
            if (limit < 0)
            {
                throw new ArgumentOutOfRangeException(nameof(limit));
            }
            int[] result = new int[limit];
            for (int i = 0; i < limit; i++)
            {
                result[i] = i;
            }
            for (int i = limit - 1; i > 0; i--)
            {
                int j = rnd.Next(i + 1);
                int swap = result[i];
                result[i] = result[j];
                result[j] = swap;
            }
            return result;
        }

        /// <summary>Formats a chromosome. Square chromosomes are shown as rows.</summary>
        public virtual string FlattenSolutiontoString(T[] sol)
        {
            if (sol == null)
            {
                return string.Empty;
            }
            StringBuilder builder = new StringBuilder();
            int perRow = (int)Math.Round(Math.Sqrt(numberofGenes));
            bool square = perRow > 0 && perRow * perRow == numberofGenes;
            for (int i = 0; i < sol.Length; i++)
            {
                if (square && i > 0 && i % perRow == 0)
                {
                    builder.Append(" | ");
                }
                builder.Append(Convert.ToString(sol[i], CultureInfo.InvariantCulture));
                builder.Append(' ');
            }
            return builder.ToString();
        }

        private void EnsureReady()
        {
            if (chromosomes == null)
            {
                throw new InvalidOperationException("Reset() must be called before running the solver.");
            }
            if (chromosomes.Length != populationSize * 3)
            {
                throw new InvalidOperationException(
                    "PopulationSize changed after Reset(). Call Reset() again before running.");
            }
            if (objectiveValue == null || objectiveValue.Length != chromosomes.Length)
            {
                throw new InvalidOperationException("The objective buffer is inconsistent. Call Reset() again.");
            }
        }

        /// <summary>Legacy optimisation enum kept for source compatibility.</summary>
        public enum Optimization
        {
            Minimization = 0,
            Maximization = 1,
            GoalMatching = 2
        }
    }
}

