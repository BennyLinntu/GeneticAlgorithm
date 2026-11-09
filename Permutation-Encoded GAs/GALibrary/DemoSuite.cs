using System;
using System.Diagnostics;
using System.Globalization;
using System.Text;

namespace GALibrary
{
    /// <summary>
    /// Console-friendly demonstration used by <c>--demo --headless</c>.
    /// It loads/uses an <see cref="CostMatrixProblem"/> and exercises the binary,
    /// permutation and real-valued encodings without creating any WinForms object.
    /// </summary>
    public static class DemoSuite
    {
        public static string RunCostMatrixDemos(CostMatrixProblem problem, int iterations, int seed, Action<string> log)
        {
            if (problem == null) throw new ArgumentNullException(nameof(problem));
            if (iterations < 1) iterations = 1;

            StringBuilder report = new StringBuilder();
            Stopwatch stopwatch = Stopwatch.StartNew();

            Write(report, log, "================================================================================");
            Write(report, log, " Genetic Algorithm Workbench - headless demo");
            Write(report, log, "================================================================================");
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Problem size : {0} x {0}", problem.Size));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Iterations   : {0}", iterations));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Random seed  : {0}", seed));
            Write(report, log, string.Empty);

            double optimum = BestPermutationByBruteForce(problem);
            if (!double.IsNaN(optimum))
            {
                Write(report, log, string.Format(CultureInfo.InvariantCulture, " Brute-force optimum : {0:0.####}", optimum));
                Write(report, log, string.Empty);
            }

            RunBinaryDemo(report, log, problem, iterations, seed);
            RunPermutationDemo(report, log, problem, iterations, seed);
            RunRealCodedDemo(report, log, iterations, seed);

            stopwatch.Stop();
            Write(report, log, string.Empty);
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Demo completed in {0:0.00} s.", stopwatch.Elapsed.TotalSeconds));
            Write(report, log, "================================================================================");
            return report.ToString();
        }

        // -----------------------------------------------------------------
        // Binary encoding
        // -----------------------------------------------------------------

        private static void RunBinaryDemo(StringBuilder report, Action<string> log, CostMatrixProblem problem, int iterations, int seed)
        {
            int size = problem.Size;
            int genes = size * size;

            Write(report, log, "--------------------------------------------------------------------------------");
            Write(report, log, " [1] Binary encoding - one-point cut, deterministic selection");
            Write(report, log, "--------------------------------------------------------------------------------");

            BinaryGA solver = new BinaryGA(
                genes,
                GAOptimizationType.Minimization,
                problem.BinaryObjective,
                BinaryCrossoverType.OnePointCut);
            solver.SetRandomSeed(seed);
            solver.PopulationSize = 30;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.05;
            solver.MutationType = GAMutationType.GeneNumberBased;
            solver.SelectionType = GASelectionType.Deterministic;
            solver.IterationLimit = iterations;
            solver.PenaltyFactor = 100.0;
            solver.LeastFinessFraction = 0.1;
            solver.Reset();
            double initial = solver.SoFartheBestObjectSolutionValue;
            solver.RunToEnd(iterations);

            WriteBinaryResult(report, log, problem, solver, initial, "one-point / deterministic");

            Write(report, log, string.Empty);
            Write(report, log, " [2] Binary encoding - N-point cut, stochastic selection");
            Write(report, log, "--------------------------------------------------------------------------------");

            BinaryGA stochastic = new BinaryGA(
                genes,
                GAOptimizationType.Minimization,
                problem.BinaryObjective,
                BinaryCrossoverType.NPointCut);
            stochastic.SetRandomSeed(seed + 1);
            stochastic.PopulationSize = 30;
            stochastic.CrossoverRate = 0.9;
            stochastic.MutationRate = 0.05;
            stochastic.MutationType = GAMutationType.ChromosomesNumberBased;
            stochastic.SelectionType = GASelectionType.Stochastic;
            stochastic.IterationLimit = iterations;
            stochastic.PenaltyFactor = 100.0;
            stochastic.LeastFinessFraction = 0.1;
            stochastic.Reset();
            double stochasticInitial = stochastic.SoFartheBestObjectSolutionValue;
            stochastic.RunToEnd(iterations);

            WriteBinaryResult(report, log, problem, stochastic, stochasticInitial, "N-point / stochastic");
            Write(report, log, string.Empty);
        }

        private static void WriteBinaryResult(
            StringBuilder report,
            Action<string> log,
            CostMatrixProblem problem,
            BinaryGA solver,
            double initial,
            string label)
        {
            int size = problem.Size;
            double raw = solver.ComputeRawObjective(solver.SoFartheBestSoltion);
            int[] violations = solver.ReturnChromosomesViolation(solver.SoFartheBestSoltion);
            int totalViolations = 0;
            for (int i = 0; i < violations.Length; i++) totalViolations += violations[i];

            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Operator set        : {0}", label));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Initial best (pen.) : {0:0.####}", initial));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Best penalised value: {0:0.####}", solver.SoFartheBestObjectSolutionValue));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Best raw setup cost : {0:0.####}", raw));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Total violations    : {0}", totalViolations));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Iterations run      : {0}", solver.IterationID));
            Write(report, log, " Best cost matrix:");
            for (int row = 0; row < size; row++)
            {
                StringBuilder line = new StringBuilder("   ");
                for (int column = 0; column < size; column++)
                {
                    line.Append(solver.SoFartheBestSoltion[row * size + column] == 0 ? "0 " : "1 ");
                }
                line.Append("  job ");
                line.Append((row + 1).ToString(CultureInfo.InvariantCulture));
                line.Append(" -> machine ");
                line.Append(AssignedMachine(solver.SoFartheBestSoltion, row, size).ToString(CultureInfo.InvariantCulture));
                Write(report, log, line.ToString());
            }
        }

        private static int AssignedMachine(byte[] genes, int job, int size)
        {
            for (int machine = 0; machine < size; machine++)
            {
                if (genes[job * size + machine] != 0) return machine + 1;
            }
            return 0;
        }

        // -----------------------------------------------------------------
        // Permutation encoding
        // -----------------------------------------------------------------

        private static void RunPermutationDemo(StringBuilder report, Action<string> log, CostMatrixProblem problem, int iterations, int seed)
        {
            Write(report, log, "--------------------------------------------------------------------------------");
            Write(report, log, " [3] Permutation encoding - PMX + inversion");
            Write(report, log, "--------------------------------------------------------------------------------");

            PermutationGA detailed = new PermutationGA(
                problem.Size,
                GAOptimizationType.Minimization,
                problem.PermutationObjective,
                PermutationCrossoverType.ParialMappedX,
                PermuatationMutationType.Inversion);
            detailed.SetRandomSeed(seed + 2);
            detailed.PopulationSize = 30;
            detailed.CrossoverRate = 0.9;
            detailed.MutationRate = 0.2;
            detailed.IterationLimit = iterations;
            detailed.Reset();
            double initial = detailed.SoFartheBestObjectSolutionValue;
            detailed.RunToEnd(iterations);

            double optimum = BestPermutationByBruteForce(problem);
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Initial best        : {0:0.####}", initial));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Best setup cost     : {0:0.####}", detailed.SoFartheBestObjectSolutionValue));
            if (!double.IsNaN(optimum))
            {
                Write(report, log, string.Format(CultureInfo.InvariantCulture, " Brute-force optimum : {0:0.####}", optimum));
                Write(report, log, string.Format(CultureInfo.InvariantCulture, " Gap to optimum      : {0:0.####}%",
                    optimum == 0 ? 0.0 : 100.0 * (detailed.SoFartheBestObjectSolutionValue - optimum) / Math.Abs(optimum)));
            }
            Write(report, log, " Best permutation    : " + FormatPermutation(detailed.SoFartheBestSoltion));
            Write(report, log, string.Empty);

            Write(report, log, " Crossover comparison (mutation = inversion):");
            foreach (PermutationCrossoverType crossover in new[]
            {
                PermutationCrossoverType.ParialMappedX,
                PermutationCrossoverType.OrderX,
                PermutationCrossoverType.PositionBasedX,
                PermutationCrossoverType.OrderBasedX,
                PermutationCrossoverType.CycleX,
                PermutationCrossoverType.SubtourExchange
            })
            {
                PermutationGA solver = BuildPermutationSolver(problem, crossover, PermuatationMutationType.Inversion, iterations, seed + 10);
                Write(report, log, string.Format(CultureInfo.InvariantCulture, "   {0,-18} best = {1:0.####}",
                    crossover, solver.SoFartheBestObjectSolutionValue));
            }

            Write(report, log, string.Empty);
            Write(report, log, " Mutation comparison (crossover = PMX):");
            foreach (PermuatationMutationType mutation in new[]
            {
                PermuatationMutationType.Inversion,
                PermuatationMutationType.Insertion,
                PermuatationMutationType.Displacement,
                PermuatationMutationType.ReciprocalExchange,
                PermuatationMutationType.Heuristic
            })
            {
                PermutationGA solver = BuildPermutationSolver(problem, PermutationCrossoverType.ParialMappedX, mutation, iterations, seed + 20);
                Write(report, log, string.Format(CultureInfo.InvariantCulture, "   {0,-18} best = {1:0.####}",
                    mutation, solver.SoFartheBestObjectSolutionValue));
            }
            Write(report, log, string.Empty);
        }

        private static PermutationGA BuildPermutationSolver(
            CostMatrixProblem problem,
            PermutationCrossoverType crossover,
            PermuatationMutationType mutation,
            int iterations,
            int seed)
        {
            PermutationGA solver = new PermutationGA(
                problem.Size,
                GAOptimizationType.Minimization,
                problem.PermutationObjective,
                crossover,
                mutation);
            solver.SetRandomSeed(seed);
            solver.PopulationSize = 30;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.2;
            solver.IterationLimit = iterations;
            solver.Reset();
            solver.RunToEnd(iterations);
            return solver;
        }

        private static string FormatPermutation(int[] permutation)
        {
            StringBuilder builder = new StringBuilder("[");
            for (int i = 0; i < permutation.Length; i++)
            {
                if (i > 0) builder.Append(", ");
                builder.Append((permutation[i] + 1).ToString(CultureInfo.InvariantCulture));
            }
            builder.Append(']');
            return builder.ToString();
        }

        // -----------------------------------------------------------------
        // Real-valued encoding
        // -----------------------------------------------------------------

        private static void RunRealCodedDemo(StringBuilder report, Action<string> log, int iterations, int seed)
        {
            Write(report, log, "--------------------------------------------------------------------------------");
            Write(report, log, " [4] Real-valued encoding - sphere function, 6 dimensions");
            Write(report, log, "--------------------------------------------------------------------------------");

            int dimensions = 6;
            double[] lower = new double[dimensions];
            double[] upper = new double[dimensions];
            for (int i = 0; i < dimensions; i++)
            {
                lower[i] = -5.12;
                upper[i] = 5.12;
            }

            RealNumberEncodedGA solver = new RealNumberEncodedGA(
                dimensions,
                lower,
                upper,
                GAOptimizationType.Minimization,
                Sphere);
            solver.SetRandomSeed(seed + 30);
            solver.PopulationSize = 40;
            solver.CrossoverRate = 0.9;
            solver.MutationRate = 0.2;
            solver.IterationLimit = Math.Max(iterations, 40);
            solver.CrossoverType = RealNumberCrossoverType.FBMS;
            solver.Reset();
            double initial = solver.SoFartheBestObjectSolutionValue;
            solver.RunToEnd(solver.IterationLimit);

            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Initial best        : {0:0.######}", initial));
            Write(report, log, string.Format(CultureInfo.InvariantCulture, " Best objective      : {0:0.######}", solver.SoFartheBestObjectSolutionValue));
            Write(report, log, " Best vector         : " + FormatDoubleVector(solver.SoFartheBestSoltion));
            Write(report, log, string.Empty);
        }

        private static double Sphere(double[] genes)
        {
            double total = 0.0;
            for (int i = 0; i < genes.Length; i++)
            {
                total += genes[i] * genes[i];
            }
            return total;
        }

        private static string FormatDoubleVector(double[] values)
        {
            StringBuilder builder = new StringBuilder("[");
            for (int i = 0; i < values.Length; i++)
            {
                if (i > 0) builder.Append(", ");
                builder.Append(values[i].ToString("0.####", CultureInfo.InvariantCulture));
            }
            builder.Append(']');
            return builder.ToString();
        }

        // -----------------------------------------------------------------
        // Exact reference for small instances
        // -----------------------------------------------------------------

        /// <summary>
        /// Exhaustive permutation search. Returns NaN for instances larger than 8
        /// (8! = 40320, still instantaneous; the cap keeps the demo predictable).
        /// </summary>
        public static double BestPermutationByBruteForce(CostMatrixProblem problem)
        {
            if (problem == null) throw new ArgumentNullException(nameof(problem));
            if (problem.Size > 8) return double.NaN;

            int[] permutation = new int[problem.Size];
            for (int i = 0; i < permutation.Length; i++) permutation[i] = i;
            double best = double.MaxValue;
            Permute(problem, permutation, 0, ref best);
            return best;
        }

        private static void Permute(CostMatrixProblem problem, int[] values, int start, ref double best)
        {
            if (start == values.Length)
            {
                double cost = problem.PermutationObjective(values);
                if (cost < best) best = cost;
                return;
            }
            for (int i = start; i < values.Length; i++)
            {
                int swap = values[start];
                values[start] = values[i];
                values[i] = swap;
                Permute(problem, values, start + 1, ref best);
                swap = values[start];
                values[start] = values[i];
                values[i] = swap;
            }
        }

        // -----------------------------------------------------------------
        // Helpers
        // -----------------------------------------------------------------

        private static void Write(StringBuilder report, Action<string> log, string message)
        {
            report.AppendLine(message);
            if (log != null) log(message);
        }
    }
}

