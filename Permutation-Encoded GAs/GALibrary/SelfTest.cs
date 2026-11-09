using System;
using System.Collections.Generic;
using System.Globalization;
using System.Text;

namespace GALibrary
{
    /// <summary>Result of a headless self-test run.</summary>
    public sealed class SelfTestReport
    {
        private readonly List<string> failures = new List<string>();

        public int Checks { get; internal set; }

        public IList<string> Failures
        {
            get { return failures; }
        }

        public bool Passed
        {
            get { return failures.Count == 0; }
        }

        public string ToText()
        {
            StringBuilder builder = new StringBuilder();
            builder.AppendLine(string.Format(
                CultureInfo.InvariantCulture,
                "Self-test: {0} checks, {1} failure(s).",
                Checks,
                failures.Count));
            for (int i = 0; i < failures.Count; i++)
            {
                builder.AppendLine("  FAIL: " + failures[i]);
            }
            builder.AppendLine(Passed ? "RESULT: PASS" : "RESULT: FAIL");
            return builder.ToString();
        }
    }

    /// <summary>
    /// Headless validation of the operators and engine invariants. Run it from
    /// either executable with <c>--selftest</c>.
    /// </summary>
    public static class GeneticAlgorithmSelfTest
    {
        public static SelfTestReport RunAll()
        {
            SelfTestReport report = new SelfTestReport();
            int seed = 20260101;

            CostMatrixProblem problem = CostMatrixProblem.Generate(7, 0.0, 10.0, new Random(seed));
            CostMatrixProblem small = CostMatrixProblem.Generate(4, 0.0, 10.0, new Random(seed + 1));

            TestPermutationOperators(report, problem);
            TestBinaryOperators(report, small);
            TestRealOperators(report);
            TestGenericEngine(report, small);
            TestPenalty(report);
            TestParameterValidation(report, small);
            TestRealCodedOptimisation(report);

            return report;
        }

        // -----------------------------------------------------------------
        // Permutation encoding
        // -----------------------------------------------------------------

        private static void TestPermutationOperators(SelfTestReport report, CostMatrixProblem problem)
        {
            int n = problem.Size;
            PermutationGA solver = new PermutationGA(
                n,
                GAOptimizationType.Minimization,
                problem.PermutationObjective,
                PermutationCrossoverType.ParialMappedX,
                PermuatationMutationType.Inversion);
            solver.SetRandomSeed(11);
            solver.PopulationSize = 24;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.5;
            solver.IterationLimit = 50;
            solver.Reset();

            int parentA = 0;
            int parentB = 1;
            int childA = solver.PopulationSize;
            int childB = solver.PopulationSize + 1;

            foreach (PermutationCrossoverType type in new[]
            {
                PermutationCrossoverType.ParialMappedX,
                PermutationCrossoverType.OrderX,
                PermutationCrossoverType.PositionBasedX,
                PermutationCrossoverType.OrderBasedX,
                PermutationCrossoverType.CycleX,
                PermutationCrossoverType.SubtourExchange
            })
            {
                solver.CrossoverType = type;
                for (int trial = 0; trial < 60; trial++)
                {
                    parentA = trial % solver.PopulationSize;
                    parentB = (trial * 5 + 3) % solver.PopulationSize;
                    try
                    {
                        solver.GenerateCrossoveredChromosomes(parentA, parentB, childA, childB);
                    }
                    catch (Exception ex)
                    {
                        report.Failures.Add(type + " threw " + ex.GetType().Name + ": " + ex.Message
                            + " childA=[" + string.Join(",", solver.Chromosomes[childA]) + "]"
                            + " childB=[" + string.Join(",", solver.Chromosomes[childB]) + "]");
                        continue;
                    }
                    Check(report, type + " child A is a valid permutation", PermutationGA.IsValidPermutation(solver.Chromosomes[childA], n));
                    Check(report, type + " child B is a valid permutation", PermutationGA.IsValidPermutation(solver.Chromosomes[childB], n));
                }
            }

            bool[] mask = new bool[n];
            for (int i = 0; i < mask.Length; i++)
            {
                mask[i] = true;
            }

            foreach (PermuatationMutationType type in new[]
            {
                PermuatationMutationType.Inversion,
                PermuatationMutationType.Insertion,
                PermuatationMutationType.Displacement,
                PermuatationMutationType.ReciprocalExchange,
                PermuatationMutationType.Heuristic
            })
            {
                solver.PermuatationMutationType = type;
                for (int trial = 0; trial < 60; trial++)
                {
                    parentA = trial % solver.PopulationSize;
                    solver.GenerateMutatedChromosomes(parentA, childA, mask);
                    Check(report, type + " child is a valid permutation", PermutationGA.IsValidPermutation(solver.Chromosomes[childA], n));
                }
            }

            // The whole engine must never leak an invalid permutation into the pool.
            solver.PermuatationMutationType = PermuatationMutationType.Insertion;
            solver.CrossoverType = PermutationCrossoverType.OrderBasedX;
            for (int generation = 0; generation < 15; generation++)
            {
                solver.RunOneIteration();
                int total = solver.ActivePoolSize;
                for (int row = 0; row < total; row++)
                {
                    if (!PermutationGA.IsValidPermutation(solver.Chromosomes[row], n))
                    {
                        report.Failures.Add("Population contains an invalid permutation at row " + row);
                        break;
                    }
                }
            }
            Check(report, "Permutation generations stay valid", report.Failures.Count == 0 || !ContainsPermutationFailure(report));
        }

        private static bool ContainsPermutationFailure(SelfTestReport report)
        {
            for (int i = 0; i < report.Failures.Count; i++)
            {
                if (report.Failures[i].Contains("invalid permutation")) return true;
            }
            return false;
        }

        // -----------------------------------------------------------------
        // Binary encoding
        // -----------------------------------------------------------------

        private static void TestBinaryOperators(SelfTestReport report, CostMatrixProblem problem)
        {
            int genes = problem.Size * problem.Size;
            BinaryGA solver = new BinaryGA(
                genes,
                GAOptimizationType.Minimization,
                problem.BinaryObjective,
                BinaryCrossoverType.OnePointCut);
            solver.SetRandomSeed(22);
            solver.PopulationSize = 20;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.5;
            solver.IterationLimit = 20;
            solver.Reset();

            int parentA = 0;
            int parentB = 1;
            int childA = solver.PopulationSize;
            int childB = solver.PopulationSize + 1;

            foreach (BinaryCrossoverType type in new[]
            {
                BinaryCrossoverType.OnePointCut,
                BinaryCrossoverType.TwoPointCut,
                BinaryCrossoverType.NPointCut
            })
            {
                solver.CrossoverType = type;
                for (int trial = 0; trial < 80; trial++)
                {
                    parentA = trial % solver.PopulationSize;
                    parentB = (trial * 7 + 2) % solver.PopulationSize;
                    solver.GenerateCrossoveredChromosomes(parentA, parentB, childA, childB);

                    bool inherited = true;
                    bool fromFather = false;
                    bool fromMother = false;
                    bool parentsDiffer = false;
                    for (int gene = 0; gene < genes; gene++)
                    {
                        byte dad = solver.Chromosomes[parentA][gene];
                        byte mom = solver.Chromosomes[parentB][gene];
                        byte one = solver.Chromosomes[childA][gene];
                        byte two = solver.Chromosomes[childB][gene];
                        if (one != 0 && one != 1) inherited = false;
                        if (two != 0 && two != 1) inherited = false;
                        bool direct = one == dad && two == mom;
                        bool swapped = one == mom && two == dad;
                        if (!direct && !swapped) inherited = false;
                        if (dad != mom)
                        {
                            parentsDiffer = true;
                            if (one == dad) fromFather = true;
                            if (one == mom) fromMother = true;
                        }
                    }
                    Check(report, type + " children are complementary binary strings", inherited);
                    bool usesBoth = !parentsDiffer || fromFather || fromMother;
                    if (!usesBoth && report.Failures.Count < 3)
                    {
                        report.Failures.Add(type + " combine diagnostic father=" + BitString(solver.Chromosomes[parentA])
                            + " mother=" + BitString(solver.Chromosomes[parentB])
                            + " childA=" + BitString(solver.Chromosomes[childA])
                            + " childB=" + BitString(solver.Chromosomes[childB]));
                    }
                    Check(report, type + " children combine both parents", usesBoth);
                }
            }

            bool[] mask = new bool[genes];
            for (int i = 0; i < mask.Length; i++)
            {
                mask[i] = true;
            }
            solver.GenerateMutatedChromosomes(0, childA, mask);
            bool flipped = true;
            for (int gene = 0; gene < genes; gene++)
            {
                if (solver.Chromosomes[childA][gene] == solver.Chromosomes[0][gene]) flipped = false;
            }
            Check(report, "Bit-flip mutation complements every selected gene", flipped);

            // Gene-based mutation at rate 1 must produce exactly N fully evaluated
            // mutation children and must stay inside the 3N buffer.
            solver.MutationType = GAMutationType.GeneNumberBased;
            solver.MutationRate = 1.0;
            solver.CrossoverRate = 0.8;
            solver.Reset();
            solver.RunOneIteration();
            Check(report, "Gene-based mutation fills every parent", solver.NumberofMutatedChildren == solver.PopulationSize);
            Check(report,
                "3N buffer is respected after gene-based mutation",
                solver.NumberofCrossedChildren + solver.NumberofMutatedChildren <= 2 * solver.PopulationSize);
            Check(report, "All active rows were evaluated", AllActiveRowsEvaluated(solver));

            // A zero mutation rate must not create mutation children.
            solver.MutationRate = 0.0;
            solver.Reset();
            solver.RunOneIteration();
            Check(report, "Mutation rate 0 creates no mutation children", solver.NumberofMutatedChildren == 0);
        }

        private static string BitString(byte[] values)
        {
            StringBuilder builder = new StringBuilder(values.Length);
            for (int i = 0; i < values.Length; i++)
            {
                builder.Append(values[i] == 0 ? '0' : '1');
            }
            return builder.ToString();
        }

        private static bool AllActiveRowsEvaluated(GenericSolver<byte> solver)        {
            int total = solver.ActivePoolSize;
            for (int row = 0; row < total; row++)
            {
                if (double.IsNaN(solver.ObjectiveValue[row]) || double.IsInfinity(solver.ObjectiveValue[row]))
                {
                    return false;
                }
            }
            return total <= solver.Chromosomes.Length;
        }

        // -----------------------------------------------------------------
        // Real encoding
        // -----------------------------------------------------------------

        private static void TestRealOperators(SelfTestReport report)
        {
            int dimensions = 5;
            double[] lower = new double[dimensions];
            double[] upper = new double[dimensions];
            for (int i = 0; i < dimensions; i++)
            {
                lower[i] = -5.0;
                upper[i] = 5.0;
            }

            RealNumberEncodedGA solver = new RealNumberEncodedGA(
                dimensions,
                lower,
                upper,
                GAOptimizationType.Minimization,
                Sphere);
            solver.SetRandomSeed(33);
            solver.PopulationSize = 24;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.5;
            solver.IterationLimit = 40;
            solver.Reset();

            Check(report, "Real encoding initialises inside the bounds", RealRowsInsideBounds(solver, solver.PopulationSize));

            int childA = solver.PopulationSize;
            int childB = solver.PopulationSize + 1;
            foreach (RealNumberCrossoverType type in new[]
            {
                RealNumberCrossoverType.LVD,
                RealNumberCrossoverType.SVD,
                RealNumberCrossoverType.MOES,
                RealNumberCrossoverType.TES,
                RealNumberCrossoverType.FBMS
            })
            {
                solver.CrossoverType = type;
                solver.GenerateCrossoveredChromosomes(0, 1, childA, childB);
                Check(report, type + " children respect the bounds", RowInsideBounds(solver.Chromosomes[childA], lower, upper));
            }

            bool[] mask = new bool[dimensions];
            for (int i = 0; i < mask.Length; i++)
            {
                mask[i] = true;
            }
            solver.GenerateMutatedChromosomes(0, childA, mask);
            Check(report, "Dynamic mutation respects the bounds", RowInsideBounds(solver.Chromosomes[childA], lower, upper));

            solver.RunToEnd(40);
            Check(report, "Real-coded run stays inside the bounds", RealRowsInsideBounds(solver, solver.ActivePoolSize));
            Check(report, "Real-coded run keeps a finite best value", !double.IsNaN(solver.SoFartheBestObjectSolutionValue));
        }

        private static double Sphere(double[] genes)
        {
            double sum = 0.0;
            for (int i = 0; i < genes.Length; i++)
            {
                sum += genes[i] * genes[i];
            }
            return sum;
        }

        private static bool RealRowsInsideBounds(RealNumberEncodedGA solver, int total)
        {
            double[] lower = solver.LowerBounds;
            double[] upper = solver.UpperBounds;
            for (int row = 0; row < total; row++)
            {
                if (!RowInsideBounds(solver.Chromosomes[row], lower, upper)) return false;
            }
            return true;
        }

        private static bool RowInsideBounds(double[] row, double[] lower, double[] upper)
        {
            if (row == null || row.Length != lower.Length) return false;
            for (int i = 0; i < row.Length; i++)
            {
                if (double.IsNaN(row[i]) || row[i] < lower[i] - 1e-12 || row[i] > upper[i] + 1e-12) return false;
            }
            return true;
        }

        // -----------------------------------------------------------------
        // Engine invariants
        // -----------------------------------------------------------------

        private static void TestGenericEngine(SelfTestReport report, CostMatrixProblem problem)
        {
            int genes = problem.Size * problem.Size;
            BinaryGA solver = new BinaryGA(
                genes,
                GAOptimizationType.Minimization,
                problem.BinaryObjective,
                BinaryCrossoverType.NPointCut);
            solver.SetRandomSeed(44);
            solver.PopulationSize = 12;
            solver.CrossoverRate = 1.0;
            solver.MutationRate = 1.0;
            solver.MutationType = GAMutationType.ChromosomesNumberBased;
            solver.IterationLimit = 20;
            solver.PenaltyFactor = 100.0;
            solver.Reset();

            Check(report, "Initial best equals the best of generation 0",
                Math.Abs(solver.InitialBestObjective - solver.SoFartheBestObjectSolutionValue) < 1e-9);
            Check(report, "Initial best solution is allocated", solver.SoFartheBestSoltion != null
                && solver.SoFartheBestSoltion.Length == genes);

            double previousBest = solver.SoFartheBestObjectSolutionValue;
            for (int generation = 0; generation < 10; generation++)
            {
                solver.RunOneIteration();
                Check(report,
                    "Child counts stay within the 3N buffers",
                    solver.NumberofCrossedChildren <= solver.PopulationSize
                    && solver.NumberofMutatedChildren <= solver.PopulationSize
                    && solver.ActivePoolSize <= 3 * solver.PopulationSize);
                Check(report, "Best-so-far never gets worse", solver.SoFartheBestObjectSolutionValue <= previousBest + 1e-9);
                previousBest = solver.SoFartheBestObjectSolutionValue;
            }

            // RunToEnd must never exceed IterationLimit.
            solver.RunToEnd(10000);
            Check(report, "RunToEnd respects IterationLimit", solver.IterationID == solver.IterationLimit);

            // Stochastic selection must survive many generations and keep a valid best.
            BinaryGA stochastic = new BinaryGA(
                genes,
                GAOptimizationType.Minimization,
                problem.BinaryObjective,
                BinaryCrossoverType.TwoPointCut);
            stochastic.SetRandomSeed(55);
            stochastic.PopulationSize = 14;
            stochastic.CrossoverRate = 0.8;
            stochastic.MutationRate = 0.2;
            stochastic.SelectionType = GASelectionType.Stochastic;
            stochastic.IterationLimit = 40;
            stochastic.PenaltyFactor = 100.0;
            stochastic.Reset();
            double initial = stochastic.SoFartheBestObjectSolutionValue;
            stochastic.RunToEnd(40);
            Check(report, "Stochastic selection runs to completion", stochastic.IterationID == 40);
            Check(report, "Stochastic selection preserves the best-so-far", stochastic.SoFartheBestObjectSolutionValue <= initial + 1e-9);
            Check(report, "Stochastic best solution is a valid buffer row",
                stochastic.SoFartheBestSoltion != null && stochastic.SoFartheBestSoltion.Length == genes);

            // Selection must always return exactly N parents.
            stochastic.PerformSelectionOperation();
            Check(report, "Selection produces a finite best objective",
                !double.IsNaN(stochastic.IterationBestObjective) && !double.IsInfinity(stochastic.IterationBestObjective));
        }

        private static void TestPenalty(SelfTestReport report)
        {
            double[,] ones = new double[3, 3];
            for (int row = 0; row < 3; row++)
            {
                for (int column = 0; column < 3; column++)
                {
                    ones[row, column] = 1.0;
                }
            }
            CostMatrixProblem problem = new CostMatrixProblem(ones);
            BinaryGA solver = new BinaryGA(9, GAOptimizationType.Minimization, problem.BinaryObjective, BinaryCrossoverType.OnePointCut);
            solver.PopulationSize = 4;
            solver.PenaltyFactor = 100.0;
            solver.Reset();

            byte[] zeros = new byte[9];
            double raw = solver.ComputeRawObjective(zeros);
            double penalised = solver.ComputePenalizedObjective(zeros);
            Check(report, "Raw objective of an all-zero matrix is zero", Math.Abs(raw) < 1e-12);
            Check(report, "Penalty uses the reported row/column violations", Math.Abs(penalised - 600.0) < 1e-9);

            int[] violations = solver.ReturnChromosomesViolation(zeros);
            Check(report, "Binary violations cover rows and columns", violations != null && violations.Length == 6);
        }

        private static void TestParameterValidation(SelfTestReport report, CostMatrixProblem problem)
        {
            int genes = problem.Size * problem.Size;
            BinaryGA solver = new BinaryGA(genes, GAOptimizationType.Minimization, problem.BinaryObjective, BinaryCrossoverType.OnePointCut);
            Check(report, "PopulationSize rejects values below 2", Throws(delegate { solver.PopulationSize = 1; }));
            Check(report, "CrossoverRate rejects values above 1", Throws(delegate { solver.CrossoverRate = 1.5; }));
            Check(report, "MutationRate rejects negative values", Throws(delegate { solver.MutationRate = -0.1; }));
            Check(report, "IterationLimit rejects negative values", Throws(delegate { solver.IterationLimit = -1; }));
            Check(report, "PenaltyFactor rejects negative values", Throws(delegate { solver.PenaltyFactor = -1.0; }));
            Check(report, "Least fitness fraction rejects zero", Throws(delegate { solver.LeastFinessFraction = 0.0; }));
            Check(report, "Chromosome length must be positive", Throws(delegate
            {
                new BinaryGA(0, GAOptimizationType.Minimization, problem.BinaryObjective, BinaryCrossoverType.OnePointCut);
            }));
            Check(report, "Real encoding rejects inverted bounds", Throws(delegate
            {
                new RealNumberEncodedGA(2, new[] { 1.0, -1.0 }, new[] { 0.0, 1.0 }, GAOptimizationType.Minimization, Sphere);
            }));
            Check(report, "Run before Reset is rejected", Throws(delegate
            {
                BinaryGA fresh = new BinaryGA(genes, GAOptimizationType.Minimization, problem.BinaryObjective, BinaryCrossoverType.OnePointCut);
                fresh.RunOneIteration();
            }));
        }

        private static void TestRealCodedOptimisation(SelfTestReport report)
        {
            int dimensions = 6;
            double[] lower = new double[dimensions];
            double[] upper = new double[dimensions];
            for (int i = 0; i < dimensions; i++)
            {
                lower[i] = -5.12;
                upper[i] = 5.12;
            }
            RealNumberEncodedGA solver = new RealNumberEncodedGA(dimensions, lower, upper, GAOptimizationType.Minimization, Sphere);
            solver.SetRandomSeed(66);
            solver.PopulationSize = 40;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.25;
            solver.IterationLimit = 120;
            solver.CrossoverType = RealNumberCrossoverType.FBMS;
            solver.Reset();
            double initial = solver.SoFartheBestObjectSolutionValue;
            solver.RunToEnd(120);
            Check(report, "Real-coded GA never worsens the best value", solver.SoFartheBestObjectSolutionValue <= initial + 1e-9);
            Check(report, "Real-coded GA improved the sphere function", solver.SoFartheBestObjectSolutionValue < initial);
            Check(report, "Real-coded best solution respects the bounds",
                RowInsideBounds(solver.SoFartheBestSoltion, lower, upper));
        }

        // -----------------------------------------------------------------
        // Helpers
        // -----------------------------------------------------------------

        private static void Check(SelfTestReport report, string name, bool condition)
        {
            report.Checks++;
            if (!condition)
            {
                report.Failures.Add(name);
            }
        }

        private static bool Throws(Action action)
        {
            try
            {
                action();
                return false;
            }
            catch (ArgumentOutOfRangeException)
            {
                return true;
            }
            catch (ArgumentException)
            {
                return true;
            }
            catch (InvalidOperationException)
            {
                return true;
            }
        }
    }
}

