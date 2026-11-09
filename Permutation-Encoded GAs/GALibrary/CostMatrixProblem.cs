using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;

namespace GALibrary
{
    /// <summary>
    /// A square selection / setup-cost matrix shared by the two WinForms
    /// front ends, the demo runner and the self-tests.
    ///
    /// File format (".aop"):
    /// <code>
    ///   n
    ///   c00 c01 ... c0(n-1)
    ///   ...
    ///   c(n-1)0 ... c(n-1)(n-1)
    /// </code>
    /// Blank lines and lines starting with '#' are ignored. Values are parsed
    /// with the invariant culture, so "1.5" works on every machine.
    /// </summary>
    public sealed class CostMatrixProblem
    {
        private readonly double[,] costs;

        public int Size { get; private set; }

        public CostMatrixProblem(double[,] costs)
        {
            if (costs == null) throw new ArgumentNullException(nameof(costs));
            int rows = costs.GetLength(0);
            int columns = costs.GetLength(1);
            if (rows != columns)
            {
                throw new ArgumentException("The cost matrix must be square.", nameof(costs));
            }
            if (rows < 2)
            {
                throw new ArgumentException("The cost matrix must be at least 2 x 2.", nameof(costs));
            }
            for (int row = 0; row < rows; row++)
            {
                for (int column = 0; column < columns; column++)
                {
                    if (double.IsNaN(costs[row, column]) || double.IsInfinity(costs[row, column]))
                    {
                        throw new ArgumentException("The cost matrix contains a non-finite value.", nameof(costs));
                    }
                }
            }

            Size = rows;
            this.costs = (double[,])costs.Clone();
        }

        public double this[int job, int machine]
        {
            get
            {
                if (job < 0 || job >= Size) throw new ArgumentOutOfRangeException(nameof(job));
                if (machine < 0 || machine >= Size) throw new ArgumentOutOfRangeException(nameof(machine));
                return costs[job, machine];
            }
        }

        public double[,] ToArray()
        {
            return (double[,])costs.Clone();
        }

        // -----------------------------------------------------------------
        // Objectives
        // -----------------------------------------------------------------

        /// <summary>Total setup cost of a row-major binary cost matrix.</summary>
        public double BinaryObjective(byte[] genes)
        {
            int expected = Size * Size;
            if (genes == null || genes.Length != expected)
            {
                throw new ArgumentException("The binary chromosome must contain n*n genes.", nameof(genes));
            }

            double total = 0.0;
            int index = 0;
            for (int job = 0; job < Size; job++)
            {
                for (int machine = 0; machine < Size; machine++)
                {
                    if (genes[index] != 0)
                    {
                        total += costs[job, machine];
                    }
                    index++;
                }
            }
            return total;
        }

        /// <summary>Total setup cost of a permutation where job i is assigned to machine perm[i].</summary>
        public double PermutationObjective(int[] permutation)
        {
            if (permutation == null || permutation.Length != Size)
            {
                throw new ArgumentException("The permutation must contain exactly n genes.", nameof(permutation));
            }

            bool[] seen = new bool[Size];
            double total = 0.0;
            for (int job = 0; job < Size; job++)
            {
                int machine = permutation[job];
                if (machine < 0 || machine >= Size || seen[machine])
                {
                    throw new ArgumentException("The chromosome is not a valid permutation.", nameof(permutation));
                }
                seen[machine] = true;
                total += costs[job, machine];
            }
            return total;
        }

        /// <summary>Row violations followed by column violations of a binary matrix.</summary>
        public int[] BinaryViolations(byte[] genes)
        {
            int expected = Size * Size;
            if (genes == null || genes.Length != expected)
            {
                throw new ArgumentException("The binary chromosome must contain n*n genes.", nameof(genes));
            }

            int[] violations = new int[2 * Size];
            for (int job = 0; job < Size; job++)
            {
                int count = 0;
                for (int machine = 0; machine < Size; machine++)
                {
                    if (genes[job * Size + machine] != 0) count++;
                }
                violations[job] = Math.Abs(count - 1);
            }
            for (int machine = 0; machine < Size; machine++)
            {
                int count = 0;
                for (int job = 0; job < Size; job++)
                {
                    if (genes[job * Size + machine] != 0) count++;
                }
                violations[Size + machine] = Math.Abs(count - 1);
            }
            return violations;
        }

        // -----------------------------------------------------------------
        // Conversions
        // -----------------------------------------------------------------

        public byte[] PermutationToBinary(int[] permutation)
        {
            PermutationObjective(permutation); // validates
            byte[] genes = new byte[Size * Size];
            for (int job = 0; job < Size; job++)
            {
                genes[job * Size + permutation[job]] = 1;
            }
            return genes;
        }

        /// <summary>
        /// Greedy decoding used for display: each job picks its cheapest still
        /// unused machine. The result is always a valid permutation but it is not
        /// guaranteed to be the selection encoded by a partial binary matrix.
        /// </summary>
        public int[] BinaryToPermutation(byte[] genes)
        {
            BinaryViolations(genes); // validates
            int[] permutation = new int[Size];
            bool[] used = new bool[Size];
            for (int job = 0; job < Size; job++)
            {
                int bestMachine = -1;
                double bestCost = double.MaxValue;
                for (int machine = 0; machine < Size; machine++)
                {
                    if (used[machine]) continue;
                    double cost = costs[job, machine];
                    if (cost < bestCost)
                    {
                        bestCost = cost;
                        bestMachine = machine;
                    }
                }
                if (bestMachine < 0)
                {
                    for (int machine = 0; machine < Size; machine++)
                    {
                        if (!used[machine]) { bestMachine = machine; break; }
                    }
                }
                permutation[job] = bestMachine;
                used[bestMachine] = true;
            }
            return permutation;
        }

        // -----------------------------------------------------------------
        // Generation and IO
        // -----------------------------------------------------------------

        public static CostMatrixProblem Generate(int size, double minimum, double maximum, Random random)
        {
            if (random == null) throw new ArgumentNullException(nameof(random));
            if (size < 2) throw new ArgumentOutOfRangeException(nameof(size), "The matrix must be at least 2 x 2.");
            if (size > 200) throw new ArgumentOutOfRangeException(nameof(size), "The matrix is limited to 200 x 200.");
            if (double.IsNaN(minimum) || double.IsNaN(maximum) || double.IsInfinity(minimum) || double.IsInfinity(maximum))
            {
                throw new ArgumentException("The minimum and maximum must be finite.");
            }
            if (minimum > maximum)
            {
                throw new ArgumentException("The minimum must not be larger than the maximum.");
            }

            double[,] data = new double[size, size];
            for (int job = 0; job < size; job++)
            {
                for (int machine = 0; machine < size; machine++)
                {
                    data[job, machine] = Math.Round(minimum + random.NextDouble() * (maximum - minimum), 3);
                }
            }
            return new CostMatrixProblem(data);
        }

        public static CostMatrixProblem Load(string path)
        {
            if (path == null) throw new ArgumentNullException(nameof(path));
            if (!File.Exists(path)) throw new FileNotFoundException("Selection file not found.", path);

            string[] rawLines = File.ReadAllLines(path);
            List<string> lines = new List<string>();
            for (int i = 0; i < rawLines.Length; i++)
            {
                string line = rawLines[i].Trim();
                if (line.Length == 0) continue;
                if (line.StartsWith("#", StringComparison.Ordinal)) continue;
                lines.Add(line);
            }

            if (lines.Count < 2)
            {
                throw new FormatException("The file must start with the matrix size followed by the cost rows.");
            }

            int size;
            if (!int.TryParse(lines[0], NumberStyles.Integer, CultureInfo.InvariantCulture, out size))
            {
                throw new FormatException("The first non-comment line must be the integer matrix size.");
            }
            if (size < 2 || size > 200)
            {
                throw new FormatException("The matrix size must be between 2 and 200.");
            }
            if (lines.Count < size + 1)
            {
                throw new FormatException(string.Format(
                    CultureInfo.InvariantCulture,
                    "Expected {0} cost rows but found {1}.",
                    size,
                    lines.Count - 1));
            }

            double[,] data = new double[size, size];
            for (int row = 0; row < size; row++)
            {
                string[] tokens = lines[row + 1].Split(new[] { ' ', '\t', ',' }, StringSplitOptions.RemoveEmptyEntries);
                if (tokens.Length < size)
                {
                    throw new FormatException(string.Format(
                        CultureInfo.InvariantCulture,
                        "Row {0} contains {1} values; {2} are required.",
                        row + 1,
                        tokens.Length,
                        size));
                }
                for (int column = 0; column < size; column++)
                {
                    double value;
                    if (!double.TryParse(tokens[column], NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out value)
                        || double.IsNaN(value) || double.IsInfinity(value))
                    {
                        throw new FormatException(string.Format(
                            CultureInfo.InvariantCulture,
                            "Row {0}, column {1} is not a finite number: '{2}'.",
                            row + 1,
                            column + 1,
                            tokens[column]));
                    }
                    data[row, column] = value;
                }
            }

            return new CostMatrixProblem(data);
        }

        public void Save(string path)
        {
            if (path == null) throw new ArgumentNullException(nameof(path));
            StringBuilder builder = new StringBuilder();
            builder.AppendLine(Size.ToString(CultureInfo.InvariantCulture));
            for (int row = 0; row < Size; row++)
            {
                for (int column = 0; column < Size; column++)
                {
                    if (column > 0) builder.Append(' ');
                    builder.Append(costs[row, column].ToString("0.###", CultureInfo.InvariantCulture));
                }
                builder.AppendLine();
            }
            File.WriteAllText(path, builder.ToString(), new UTF8Encoding(false));
        }
    }
}

