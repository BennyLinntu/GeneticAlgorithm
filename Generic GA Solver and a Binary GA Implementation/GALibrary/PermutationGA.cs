using System;
using System.Collections.Generic;
using System.Linq;

namespace GALibrary
{
    /// <summary>Crossover operators for permutation-encoded chromosomes.</summary>
    public enum PermutationCrossoverType
    {
        /// <summary>Partially matched crossover (PMX).</summary>
        ParialMappedX,

        /// <summary>Order crossover (OX).</summary>
        OrderX,

        /// <summary>Position-based crossover (POS).</summary>
        PositionBasedX,

        /// <summary>Order-based crossover (OX2).</summary>
        OrderBasedX,

        /// <summary>Cycle crossover (CX).</summary>
        CycleX,

        /// <summary>Subtour exchange with duplicate repair.</summary>
        SubtourExchange
    }

    /// <summary>Mutation operators for permutation-encoded chromosomes.</summary>
    public enum PermuatationMutationType
    {
        Inversion,
        Insertion,
        Displacement,
        ReciprocalExchange,
        Heuristic
    }

    /// <summary>
    /// Permutation encoding: every chromosome is a permutation of 0..n-1, which
    /// makes it a natural representation for selection, sequencing and routing
    /// problems. All operators below are permutation-preserving.
    /// </summary>
    public class PermutationGA : GenericSolver<int>
    {
        public PermutationCrossoverType CrossoverType { get; set; } = PermutationCrossoverType.ParialMappedX;

        public PermuatationMutationType PermuatationMutationType { get; set; } = PermuatationMutationType.Inversion;

        public PermutationGA(
            int numberofGenes,
            GAOptimizationType optimizationType,
            ObjectiveFunction<int> objectiveFunction,
            PermutationCrossoverType crossoverType,
            PermuatationMutationType mutationType)
            : base(numberofGenes, optimizationType, objectiveFunction)
        {
            if (numberofGenes < 2)
            {
                throw new ArgumentOutOfRangeException(nameof(numberofGenes), "A permutation requires at least two genes.");
            }
            CrossoverType = crossoverType;
            PermuatationMutationType = mutationType;
        }

        // -----------------------------------------------------------------
        // Initialisation
        // -----------------------------------------------------------------

        public override void InitializePopulation()
        {
            for (int row = 0; row < PopulationSize; row++)
            {
                int[] permutation = CreateRandomPermutation(numberofGenes);
                Array.Copy(permutation, Chromosomes[row], numberofGenes);
                ObjectiveValue[row] = Evaluate(Chromosomes[row]);
            }
        }

        /// <summary>Returns a shuffled array of <paramref name="n"/> values, taking the first <paramref name="take"/>.</summary>
        public int[] ShuffleaArray(int n, int take)
        {
            if (n < 0) throw new ArgumentOutOfRangeException(nameof(n));
            if (take < 0 || take > n) throw new ArgumentOutOfRangeException(nameof(take));
            int[] values = ShuffleIndiceArray(n);
            if (take == n)
            {
                return values;
            }
            int[] result = new int[take];
            Array.Copy(values, result, take);
            return result;
        }

        private int[] CreateRandomPermutation(int length)
        {
            return ShuffleIndiceArray(length);
        }

        /// <summary>True when <paramref name="chromosome"/> is a permutation of 0..n-1.</summary>
        public static bool IsValidPermutation(int[] chromosome, int length)
        {
            if (chromosome == null || chromosome.Length != length) return false;
            bool[] seen = new bool[length];
            for (int i = 0; i < length; i++)
            {
                int value = chromosome[i];
                if (value < 0 || value >= length || seen[value]) return false;
                seen[value] = true;
            }
            return true;
        }

        private int[] RandomCutPair()
        {
            int first = rnd.Next(0, numberofGenes - 1);
            int second = rnd.Next(first + 1, numberofGenes);
            return new[] { first, second };
        }

        private int[] RandomPositions(int count)
        {
            if (count < 0 || count > numberofGenes) throw new ArgumentOutOfRangeException(nameof(count));
            int[] all = ShuffleIndiceArray(numberofGenes);
            int[] result = new int[count];
            Array.Copy(all, result, count);
            return result;
        }

        // -----------------------------------------------------------------
        // Crossover operators
        // -----------------------------------------------------------------

        /// <summary>
        /// Partially matched crossover. The child copies a contiguous segment
        /// from one parent; the remaining values are taken from the other parent
        /// and translated through the position-wise segment mapping.
        /// </summary>
        private void PartialMappedCrossover(int father, int mother, int childA, int childB)
        {
            int length = numberofGenes;
            int[] cuts = RandomCutPair();
            int start = cuts[0];
            int end = cuts[1];

            int[] dad = Chromosomes[father];
            int[] mom = Chromosomes[mother];
            int[] childOne = new int[length];
            int[] childTwo = new int[length];
            for (int i = 0; i < length; i++)
            {
                childOne[i] = -1;
                childTwo[i] = -1;
            }

            bool[] inMomSegment = new bool[length];
            bool[] inDadSegment = new bool[length];
            int[] momToDad = new int[length];
            int[] dadToMom = new int[length];
            for (int i = 0; i < length; i++)
            {
                momToDad[i] = -1;
                dadToMom[i] = -1;
            }

            for (int i = start; i < end; i++)
            {
                childOne[i] = mom[i];
                childTwo[i] = dad[i];
                inMomSegment[mom[i]] = true;
                inDadSegment[dad[i]] = true;
                momToDad[mom[i]] = dad[i];
                dadToMom[dad[i]] = mom[i];
            }

            for (int i = 0; i < length; i++)
            {
                if (i >= start && i < end) continue;

                int valueOne = dad[i];
                int guardOne = 0;
                while (inMomSegment[valueOne] && guardOne++ <= length)
                {
                    valueOne = momToDad[valueOne];
                }
                childOne[i] = valueOne;

                int valueTwo = mom[i];
                int guardTwo = 0;
                while (inDadSegment[valueTwo] && guardTwo++ <= length)
                {
                    valueTwo = dadToMom[valueTwo];
                }
                childTwo[i] = valueTwo;
            }

            Array.Copy(childOne, Chromosomes[childA], length);
            Array.Copy(childTwo, Chromosomes[childB], length);
        }

        /// <summary>
        /// Order crossover. The child inherits a contiguous segment from one
        /// parent and fills the remaining positions in the order in which they
        /// appear in the other parent.
        /// </summary>
        private void OrderCrossover(int father, int mother, int childA, int childB)
        {
            int length = numberofGenes;
            int[] cuts = RandomCutPair();
            int start = cuts[0];
            int end = cuts[1];

            int[] childOne = OrderChild(Chromosomes[mother], Chromosomes[father], start, end);
            int[] childTwo = OrderChild(Chromosomes[father], Chromosomes[mother], start, end);

            Array.Copy(childOne, Chromosomes[childA], length);
            Array.Copy(childTwo, Chromosomes[childB], length);
        }

        private int[] OrderChild(int[] segmentParent, int[] orderParent, int start, int end)
        {
            int length = numberofGenes;
            int[] child = new int[length];
            bool[] used = new bool[length];

            for (int i = start; i < end; i++)
            {
                child[i] = segmentParent[i];
                used[child[i]] = true;
            }

            int fill = end % length;
            for (int offset = 0; offset < length; offset++)
            {
                int source = (end + offset) % length;
                int value = orderParent[source];
                if (used[value]) continue;
                child[fill] = value;
                used[value] = true;
                fill = (fill + 1) % length;
            }
            return child;
        }

        /// <summary>
        /// Position-based crossover. A random set of positions is copied from one
        /// parent and the remaining positions are filled with the other parent's
        /// values in their original order.
        /// </summary>
        private void PositionBasedCrossover(int father, int mother, int childA, int childB)
        {
            int length = numberofGenes;
            int count = rnd.Next(1, length);
            int[] positions = RandomPositions(count);
            bool[] selected = new bool[length];
            for (int i = 0; i < positions.Length; i++)
            {
                selected[positions[i]] = true;
            }

            int[] childOne = PositionChild(Chromosomes[father], Chromosomes[mother], selected);
            int[] childTwo = PositionChild(Chromosomes[mother], Chromosomes[father], selected);

            Array.Copy(childOne, Chromosomes[childA], length);
            Array.Copy(childTwo, Chromosomes[childB], length);
        }

        private int[] PositionChild(int[] selectedParent, int[] fillerParent, bool[] selected)
        {
            int length = numberofGenes;
            int[] child = new int[length];
            bool[] used = new bool[length];
            for (int i = 0; i < length; i++)
            {
                if (selected[i])
                {
                    child[i] = selectedParent[i];
                    used[child[i]] = true;
                }
            }
            int fill = 0;
            for (int i = 0; i < length; i++)
            {
                int value = fillerParent[i];
                if (used[value]) continue;
                while (fill < length && selected[fill]) fill++;
                if (fill >= length) break;
                child[fill] = value;
                used[value] = true;
                fill++;
            }
            return child;
        }

        /// <summary>
        /// Order-based crossover. A random set of positions is taken from one
        /// parent and the other parent's remaining values are kept in relative order.
        /// </summary>
        private void OrderBasedCrossover(int father, int mother, int childA, int childB)
        {
            int length = numberofGenes;
            int count = rnd.Next(1, length);
            int[] positions = RandomPositions(count);
            bool[] selected = new bool[length];
            for (int i = 0; i < positions.Length; i++)
            {
                selected[positions[i]] = true;
            }

            int[] childOne = OrderBasedChild(Chromosomes[father], Chromosomes[mother], selected);
            int[] childTwo = OrderBasedChild(Chromosomes[mother], Chromosomes[father], selected);

            Array.Copy(childOne, Chromosomes[childA], length);
            Array.Copy(childTwo, Chromosomes[childB], length);
        }

        private int[] OrderBasedChild(int[] selectedParent, int[] fillerParent, bool[] selected)
        {
            int length = numberofGenes;
            int[] child = new int[length];
            bool[] used = new bool[length];
            for (int i = 0; i < length; i++)
            {
                if (!selected[i]) continue;
                child[i] = selectedParent[i];
                used[child[i]] = true;
            }
            int fill = 0;
            for (int i = 0; i < length; i++)
            {
                int value = fillerParent[i];
                if (used[value]) continue;
                while (fill < length && selected[fill]) fill++;
                if (fill >= length) break;
                child[fill] = value;
                used[value] = true;
                fill++;
            }
            return child;
        }

        /// <summary>
        /// Cycle crossover. The father-to-mother position mapping is decomposed
        /// into cycles; each child copies whole cycles from alternating parents.
        /// </summary>
        private void CycleCrossover(int father, int mother, int childA, int childB)
        {
            int length = numberofGenes;
            int[] dad = Chromosomes[father];
            int[] mom = Chromosomes[mother];

            int[] positionInMom = new int[length];
            for (int i = 0; i < length; i++)
            {
                positionInMom[mom[i]] = i;
            }

            int[] childOne = new int[length];
            int[] childTwo = new int[length];
            bool[] visited = new bool[length];

            int cycleIndex = 0;
            for (int start = 0; start < length; start++)
            {
                if (visited[start]) continue;

                List<int> cycle = new List<int>();
                int current = start;
                while (!visited[current])
                {
                    visited[current] = true;
                    cycle.Add(current);
                    current = positionInMom[dad[current]];
                }

                bool fromFather = (cycleIndex % 2) == 0;
                for (int i = 0; i < cycle.Count; i++)
                {
                    int position = cycle[i];
                    if (fromFather)
                    {
                        childOne[position] = dad[position];
                        childTwo[position] = mom[position];
                    }
                    else
                    {
                        childOne[position] = mom[position];
                        childTwo[position] = dad[position];
                    }
                }
                cycleIndex++;
            }

            Array.Copy(childOne, Chromosomes[childA], length);
            Array.Copy(childTwo, Chromosomes[childB], length);
        }

        /// <summary>
        /// Subtour exchange. A random window is copied from each parent into the
        /// other parent's child, then duplicate values are repaired with the
        /// missing values so that both children remain valid permutations.
        /// </summary>
        private void SubtourExchangeCrossover(int father, int mother, int childA, int childB)
        {
            int length = numberofGenes;
            int windowLength = rnd.Next(1, Math.Max(2, length / 2 + 1));
            if (windowLength >= length) windowLength = length - 1;
            int startOne = rnd.Next(0, length - windowLength + 1);
            int startTwo = rnd.Next(0, length - windowLength + 1);

            int[] childOne = new int[length];
            int[] childTwo = new int[length];
            Array.Copy(Chromosomes[father], childOne, length);
            Array.Copy(Chromosomes[mother], childTwo, length);

            for (int i = 0; i < windowLength; i++)
            {
                childOne[startTwo + i] = Chromosomes[mother][startTwo + i];
                childTwo[startOne + i] = Chromosomes[father][startOne + i];
            }

            RepairDuplicates(childOne, startTwo, windowLength);
            RepairDuplicates(childTwo, startOne, windowLength);

            Array.Copy(childOne, Chromosomes[childA], length);
            Array.Copy(childTwo, Chromosomes[childB], length);
        }

        private void RepairDuplicates(int[] child, int protectedStart, int protectedLength)
        {
            int length = child.Length;
            int[] counts = new int[length];
            for (int i = 0; i < length; i++)
            {
                if (child[i] >= 0 && child[i] < length) counts[child[i]]++;
            }

            List<int> missing = new List<int>();
            for (int value = 0; value < length; value++)
            {
                if (counts[value] == 0) missing.Add(value);
            }
            int missingIndex = 0;

            bool[] kept = new bool[length];
            int protectedEnd = protectedStart + protectedLength;

            // Reserve the copied window first so an earlier occurrence of the
            // same value outside the window is recognised as a duplicate.
            for (int i = protectedStart; i < protectedEnd; i++)
            {
                kept[child[i]] = true;
            }

            for (int i = 0; i < length; i++)
            {
                if (i >= protectedStart && i < protectedEnd) continue;
                int value = child[i];
                if (value < 0 || value >= length || kept[value])
                {
                    if (missingIndex < missing.Count)
                    {
                        child[i] = missing[missingIndex];
                        kept[child[i]] = true;
                        missingIndex++;
                    }
                }
                else
                {
                    kept[value] = true;
                }
            }
        }

        // -----------------------------------------------------------------
        // Mutation operators
        // -----------------------------------------------------------------

        private void InversionMutation(int beforeMutation, int afterMutation, bool[] mutatedFlag)
        {
            int[] cuts = RandomCutPair();
            int start = cuts[0];
            int end = cuts[1];
            int[] source = Chromosomes[beforeMutation];
            int[] destination = Chromosomes[afterMutation];

            for (int i = 0; i < numberofGenes; i++)
            {
                if (i >= start && i < end)
                {
                    destination[i] = source[end - 1 - (i - start)];
                }
                else
                {
                    destination[i] = source[i];
                }
            }
        }

        private void InsertionMutation(int beforeMutation, int afterMutation, bool[] mutatedFlag)
        {
            int from = rnd.Next(numberofGenes);
            int to = rnd.Next(numberofGenes);
            InsertShift(Chromosomes[beforeMutation], Chromosomes[afterMutation], from, to);
        }

        private void DisplacementMutation(int beforeMutation, int afterMutation, bool[] mutatedFlag)
        {
            int length = numberofGenes;
            int windowLength = rnd.Next(1, length);
            int start = rnd.Next(0, length - windowLength + 1);
            int insertAt = rnd.Next(0, length - windowLength + 1);

            int[] source = Chromosomes[beforeMutation];
            int[] destination = Chromosomes[afterMutation];

            int[] remaining = new int[length - windowLength];
            int r = 0;
            for (int i = 0; i < length; i++)
            {
                if (i >= start && i < start + windowLength) continue;
                remaining[r] = source[i];
                r++;
            }

            int write = 0;
            for (int i = 0; i < remaining.Length; i++)
            {
                if (write == insertAt)
                {
                    for (int k = 0; k < windowLength; k++)
                    {
                        destination[write] = source[start + k];
                        write++;
                    }
                }
                destination[write] = remaining[i];
                write++;
            }
            if (write == insertAt)
            {
                for (int k = 0; k < windowLength; k++)
                {
                    destination[write] = source[start + k];
                    write++;
                }
            }
        }

        private void ReciprocalExchangeMutation(int beforeMutation, int afterMutation, bool[] mutatedFlag)
        {
            int first = rnd.Next(numberofGenes);
            int second = rnd.Next(numberofGenes - 1);
            if (second >= first) second++;
            Array.Copy(Chromosomes[beforeMutation], Chromosomes[afterMutation], numberofGenes);
            int swap = Chromosomes[afterMutation][first];
            Chromosomes[afterMutation][first] = Chromosomes[afterMutation][second];
            Chromosomes[afterMutation][second] = swap;
        }

        /// <summary>
        /// Heuristic mutation: hill-climb over a handful of random insertion
        /// moves and keep the best permutation found. The child always remains a
        /// valid permutation, even when no improving move exists.
        /// </summary>
        private void HeuristicMutation(int beforeMutation, int afterMutation, bool[] mutatedFlag)
        {
            int length = numberofGenes;
            int[] best = new int[length];
            Array.Copy(Chromosomes[beforeMutation], best, length);
            double bestObjective = Evaluate(best);

            int trials = Math.Max(4, Math.Min(32, length * 2));
            for (int trial = 0; trial < trials; trial++)
            {
                int from = rnd.Next(length);
                int to = rnd.Next(length);
                if (from == to) continue;

                int[] candidate = new int[length];
                InsertShift(best, candidate, from, to);
                double candidateObjective = Evaluate(candidate);
                bool better = OptimizationType == GAOptimizationType.Maximization
                    ? candidateObjective > bestObjective
                    : candidateObjective < bestObjective;
                if (better)
                {
                    Array.Copy(candidate, best, length);
                    bestObjective = candidateObjective;
                }
            }

            Array.Copy(best, Chromosomes[afterMutation], length);
        }

        private static void InsertShift(int[] source, int[] destination, int from, int to)
        {
            int length = source.Length;
            if (from == to)
            {
                Array.Copy(source, destination, length);
                return;
            }
            int[] buffer = new int[length];
            int index = 0;
            for (int i = 0; i < length; i++)
            {
                if (i == from) continue;
                buffer[index] = source[i];
                index++;
            }
            index = 0;
            for (int i = 0; i < length; i++)
            {
                if (i == to)
                {
                    destination[i] = source[from];
                }
                else
                {
                    destination[i] = buffer[index];
                    index++;
                }
            }
        }

        // -----------------------------------------------------------------
        // Dispatch
        // -----------------------------------------------------------------

        public override void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb)
        {
            ValidateIndices(father, mother, childa, childb);
            switch (CrossoverType)
            {
                case PermutationCrossoverType.ParialMappedX:
                    PartialMappedCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.OrderX:
                    OrderCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.PositionBasedX:
                    PositionBasedCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.OrderBasedX:
                    OrderBasedCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.CycleX:
                    CycleCrossover(father, mother, childa, childb);
                    break;
                case PermutationCrossoverType.SubtourExchange:
                    SubtourExchangeCrossover(father, mother, childa, childb);
                    break;
                default:
                    throw new ArgumentOutOfRangeException("CrossoverType", CrossoverType, "Unknown crossover operator.");
            }

            ObjectiveValue[childa] = Evaluate(Chromosomes[childa]);
            ObjectiveValue[childb] = Evaluate(Chromosomes[childb]);
        }

        public override void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            ValidateIndices(beforemutation, aftermutation);
            switch (PermuatationMutationType)
            {
                case PermuatationMutationType.Inversion:
                    InversionMutation(beforemutation, aftermutation, mutatedFlag);
                    break;
                case PermuatationMutationType.Insertion:
                    InsertionMutation(beforemutation, aftermutation, mutatedFlag);
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
                    throw new ArgumentOutOfRangeException("PermuatationMutationType", PermuatationMutationType, "Unknown mutation operator.");
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

