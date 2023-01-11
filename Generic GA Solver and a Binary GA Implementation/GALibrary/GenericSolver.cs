using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;

namespace GALibrary
{
    // here is enum
    public enum GAOptimizationType { Minimization, Maximization }
    public enum GAMutationType { GeneNumberBased, ChromosomesNumberBased }
    public enum GASelectionType { Deterministic, Stochastic }

    public delegate double COPObjectiveFunction(double[] solution);
    public class GenericSolver<T>
    {
        // here we show the data fileds
        public Random rnd = new Random();

        protected T[][] chromosomes;
        protected double[] objectiveValue;
        protected double[] fitnessValue;

        T[] soFartheBestSoltion;
        double soFartheBestObjectSolutionValue;
        double initialBestObjective;
        double iterationBestObjective;
        double iterationAverageObjective;

        protected int[] indices;
        protected T[][] selectedChromosome;
        protected double[] selectedObjectiveValue;
        protected int numberofGenes;

        protected ObjectiveFunction<T> objectivefunction = null;
        protected double lessFitnessFraction;

        // initial data
        private int populationSize = 10;
        private double crossoverRate = 0.8;
        private double mutationRate = 0.2;
        double leastFinessFraction = 0.1;
        double penaltyFactor = 100;

        GAOptimizationType optimizationType = GAOptimizationType.Minimization;
        GAMutationType mutationType = GAMutationType.ChromosomesNumberBased;
        GASelectionType selectionType = GASelectionType.Deterministic;

        public delegate double ObjectiveFunction<S>(S[] Solution);

        int iterationLimit = 200;
        int numberofCrossedChildren;
        int numberofMutatedChildren;

        Series seriesSoFartheBest;
        Series seriesAverage;
        Series seriesIterationtheBest;
        Optimization op = Optimization.Minimization;
        COPObjectiveFunction thobjectiveFunction;
        
        bool[][] mutatedGenes;

        int iterationID;
        // Property


        public int PopulationSize
        { 
            get => populationSize;
            set
            {
                if(value > 1)
                {
                    populationSize = value;
                }
            }
        }
        public double CrossoverRate 
        {
            get => crossoverRate;
            set
            {
                if(value <= 1 && value >= 0)
                {
                    crossoverRate = value;
                }
            }
        }
        public double MutationRate { get => mutationRate; set => mutationRate = value; }
        public T[][] Chromosomes { get => chromosomes; set => chromosomes = value; }
        public GASelectionType SelectionType { get => selectionType; set => selectionType = value; }
        public T[] SoFartheBestSoltion { get => soFartheBestSoltion; set => soFartheBestSoltion = value; }
        public double SoFartheBestObjectSolutionValue { get => soFartheBestObjectSolutionValue; set => soFartheBestObjectSolutionValue = value; }
        public int IterationID { get => iterationID; set => iterationID = value; }
        public int NumberofCrossedChildren { get => numberofCrossedChildren; set => numberofCrossedChildren = value; }
        public int NumberofMutatedChildren { get => numberofMutatedChildren; set => numberofMutatedChildren = value; }
        public Series SeriesSoFartheBest { get => seriesSoFartheBest; set => seriesSoFartheBest = value; }
        public Series SeriesAverage { get => seriesAverage; set => seriesAverage = value; }
        public Series SeriesIterationtheBest { get => seriesIterationtheBest; set => seriesIterationtheBest = value; }
        public int IterationLimit { get => iterationLimit; set => iterationLimit = value; }
        public double PenaltyFactor { get => penaltyFactor; set => penaltyFactor = value; }
        public double LeastFinessFraction { get => leastFinessFraction; set => leastFinessFraction = value; }
        public GAOptimizationType OptimizationType { get => optimizationType; set => optimizationType = value; }
        public GAMutationType MutationType { get => mutationType; set => mutationType = value; }
        //public double[] ObjectiveValue { get => ObjectiveValue1; set => ObjectiveValue1 = value; }
        public double[] ObjectiveValue1 { get => objectiveValue; set => objectiveValue = value; }
        public double[] ObjectiveValue { get => objectiveValue; set => objectiveValue = value; }
        public double InitialBestObjective { get => initialBestObjective; set => initialBestObjective = value; }



        // function

        public GenericSolver(int numberofGenes, GAOptimizationType optimizationType, ObjectiveFunction<T> objectiveFunction)
        {
            this.numberofGenes = numberofGenes;
            this.optimizationType = optimizationType;
            this.objectivefunction = objectiveFunction;
            

            // set the series

            SeriesSoFartheBest = new Series("So Far the Best Solution");
            SeriesSoFartheBest.ChartType = SeriesChartType.Line;
            SeriesSoFartheBest.Color = Color.Red;
            SeriesSoFartheBest.BorderWidth = 2;

            SeriesAverage = new Series("Solution Average");
            SeriesAverage.ChartType = SeriesChartType.Line;
            SeriesAverage.Color = Color.Green;
            SeriesAverage.BorderWidth = 2;

            SeriesIterationtheBest = new Series("Iteration the Best Soltion");
            SeriesIterationtheBest.ChartType = SeriesChartType.Line;
            SeriesSoFartheBest.Color = Color.Blue;
            SeriesIterationtheBest.BorderWidth = 2;
        }


        public void SetFitnessandObjectives(int total, double leastFitnessFraction)
        {

            double omin = double.MaxValue;
            double omax = double.MinValue;

            for (int i = 0; i < total; i++)
            {
                if(omax <objectiveValue[i])
                {
                    omax = objectiveValue[i];
                }
                if(omin > objectiveValue[i])
                {
                    omin = objectiveValue[i];
                }
            }
            double beta = Math.Max(leastFitnessFraction * (omax - omin), 1e-5);

            switch(OptimizationType)
            {
                case GAOptimizationType.Maximization:
                    for (int i = 0; i < total; i++)
                    {
                        fitnessValue[i] = beta + (objectiveValue[i] - omin);
                    }
                    break;
                case GAOptimizationType.Minimization:
                    for (int i = 0; i < total; i++)
                    {
                        fitnessValue[i] = beta + (omax - objectiveValue[i]);
                    }
                    break;
                default:
                    break;
            }
        }
   

        void CopyAllSelectionChromosomes()
        {

            for (int i = 0; i < populationSize; i++)
            {
                for (int j = 0; j < numberofGenes; j++)
                {
                    chromosomes[i][j] = selectedChromosome[i][j];
                }
                objectiveValue[i] = selectedObjectiveValue[i];
            }
        }

        void CopyThisChromestoaSelection(int Chroindex, int Seleindex)
        {
            for (int j = 0; j < numberofGenes; j++)
            {
                selectedChromosome[Seleindex][j] = chromosomes[Chroindex][j];
            }
            selectedObjectiveValue[Seleindex] = objectivefunction(selectedChromosome[Seleindex]);
        }

        public void PerformSelectionOperation()
        {
            int total = populationSize + NumberofCrossedChildren + NumberofMutatedChildren;
            SetFitnessandObjectives(total, LeastFinessFraction);

            indices = new int[total];
            for (int i = 0; i < total; i++)
            {
                indices[i] = i;
            }
            Array.Sort(fitnessValue, indices, 0, total);
            Array.Reverse(fitnessValue, 0, total);
            Array.Reverse(indices, 0, total);

            if(SelectionType == GASelectionType.Stochastic)
            {
                double Ratio = 1.0 / fitnessValue.Sum();
                double[] dartboard = new double[total];
                int[] location = new int[populationSize];
                double sum = 0;
                int dartsCounts = 0;
                double dart = rnd.NextDouble();

                for (int i = 0; i < total; i++)
                {
                    sum += Ratio * fitnessValue[i] * populationSize;
                    dartboard[i] = sum;
                }
                for (int p = 0; p < total; p++)
                {
                    while(dartboard[p] > dart && dartsCounts < populationSize)
                    {
                        location[dartsCounts] = p;
                        dartsCounts++;
                        dart++;
                    }
                }

                for (int l = 0; l < location.Length; l++)
                {
                    CopyThisChromestoaSelection(indices[location[l]], l);
                }
            }
            else
            {
                for (int i = 0; i < populationSize; i++)
                {
                    CopyThisChromestoaSelection(indices[i], i);
                }
            }

            //FittoSelectionThenDropUnfit();
        }

        public virtual int[] ReturnChromosomesViolation(T[] chrom)
        {
            throw new NotImplementedException();
        }
        private void UpdateSoFartheBestandIterationtheBest()
        {
            int total = populationSize + numberofCrossedChildren + numberofMutatedChildren;
            int the_Best_Index = 0;
            if (optimizationType == GAOptimizationType.Maximization)
            {
                for (int i = 0; i < total; i++)
                {
                    if (objectiveValue[i] > soFartheBestObjectSolutionValue)
                    {
                        // change so far the best 
                        soFartheBestObjectSolutionValue = ObjectiveValue1[i];
                        the_Best_Index = i;
                    }
                }
            }
            else
            {
                for (int i = 0; i < total; i++)
                {
                    if (objectiveValue[i] < soFartheBestObjectSolutionValue)
                    {
                        // change so far the best 
                        soFartheBestObjectSolutionValue = ObjectiveValue1[i];
                        the_Best_Index = i;
                    }
                }
            }

            // copy the solution
            for (int i = 0; i < numberofGenes; i++)
            {
                soFartheBestSoltion[i] = chromosomes[the_Best_Index][i];
            }

            // update iteration the best
            if (optimizationType == GAOptimizationType.Maximization)
            {
                iterationBestObjective = ObjectiveValue1[the_Best_Index];
            }
                
            else
            {
                iterationBestObjective = ObjectiveValue1[the_Best_Index];
            }
                
            iterationAverageObjective = (ObjectiveValue1.Sum() / total);

        }

        //public void FittoSelectionThenDropUnfit()
        //{
        //    double selectedMin = selectedObjectiveValue.Min();
        //    double selectedMax = selectedObjectiveValue.Max();
        //    int total = populationSize + NumberofCrossedChildren + NumberofMutatedChildren;
        //    switch (OptimizationType)
        //    {
        //        case GAOptimizationType.Minimization:
        //            if (SoFartheBestObjectSolutionValue > selectedMin)
        //            {
        //                SoFartheBestObjectSolutionValue = selectedMin;
        //                SoFartheBestSoltion = selectedChromosome[0];
        //            }
        //            break;
        //        case GAOptimizationType.Maximization:
        //            if (SoFartheBestObjectSolutionValue < selectedMax)
        //            {
        //                SoFartheBestObjectSolutionValue = selectedMax;
        //                SoFartheBestSoltion = selectedChromosome[0];
        //            }
        //            break;
        //        default:
        //            break;
        //    }
        //    CopyAllSelectionChromosomes();
        //    if(OptimizationType == GAOptimizationType.Maximization)
        //    {
        //        iterationAverageObjective = selectedMax;
        //    }
        //    else
        //    {
        //        iterationAverageObjective = selectedMin;
        //    }
        //    iterationAverageObjective = (objectiveValue.Sum() / total);

        //}


        public virtual void InitializePopulation()
        {
            throw new NotImplementedException();
        }
        public virtual void ClearChromsomesandObjectiveArray()
        {
            throw new Exception("");
        }

        protected void PerformCrossoverOperation()
        {
            indices = ShuffleIndiceArray(PopulationSize);
            NumberofCrossedChildren = (int)(populationSize * CrossoverRate) / 2 * 2;
            if (NumberofCrossedChildren % 2 == 1) NumberofCrossedChildren--;
            int pair = NumberofCrossedChildren / 2;

            int father, mother, childa, childb;
            for (int p = 0; p < pair; p++)
            {
                int crossoverIDa = populationSize + 2 * p;
                int crossoverIDb = populationSize + 2 * p + 1;
                father = indices[p * 2];
                mother = indices[p * 2 + 1];
                childa = crossoverIDa;
                childb = crossoverIDb;

                GenerateCrossoveredChromosomes(father, mother, childa, childb);
            }
        }


        protected void PerformMutationOperation()
        {
            for (int i = 0; i < populationSize; i++)
            {
                indices[i] = i;
                for (int g = 0; g < numberofGenes; g++)
                {
                    mutatedGenes[i][g] = false;
                }
            }
            if (MutationType == GAMutationType.GeneNumberBased)
            {
                int genespool = populationSize * numberofGenes;
                int numberofMutatedGenesPools = (int)(MutationRate * genespool);

                for (int i = 0; i < numberofMutatedGenesPools; i++)
                {
                    int genePoolLocation = rnd.Next(genespool);
                    int rowID, colID;
                    rowID = genePoolLocation / numberofGenes;
                    colID = genePoolLocation % numberofGenes;
                    indices[rowID] = int.MinValue;
                    mutatedGenes[rowID][colID] = true;
                }

                NumberofMutatedChildren = 0;
                for (int i = 0; i < populationSize; i++)
                {
                    if (indices[i] == int.MinValue)
                    {
                        GenerateMutatedChromosomes(i, populationSize + NumberofCrossedChildren + i, mutatedGenes[i]);
                    }
                    NumberofMutatedChildren++;
                }
            }
            else
            {
                indices = ShuffleIndiceArray(populationSize);
                NumberofMutatedChildren = (int)(populationSize * MutationRate);
                for (int i = 0; i < NumberofMutatedChildren; i++)
                {
                    int mutatedID = populationSize + NumberofCrossedChildren + i;
                    int geneposition = rnd.Next(numberofGenes);
                    mutatedGenes[indices[i]][geneposition] = true;

                    GenerateMutatedChromosomes(indices[i], mutatedID, mutatedGenes[indices[i]]);
                }
            }
        }

        public int[] ShuffleIndiceArray(int limit)
        {
            int[] arr = new int[limit];
            Random r = new Random();
            for (int i = 0; i < limit; i++)
            {
                arr[i] = i;
            }
            arr = arr.OrderBy(x => r.Next()).ToArray();
            return arr;
        }
      


        public virtual void GenerateCrossoveredChromosomes(int father, int mother, int childa, int childb)
        {
            throw new NotImplementedException();
        }
        public virtual void GenerateMutatedChromosomes(int beforemutation, int aftermutation, bool[] mutatedFlag)
        {
            throw new NotImplementedException();
        }


        public virtual string FlattenSolutiontoString(T[] sol)
        {
            string str = " ";
            double rowEnd;
            for (int i = 0; i < sol.Length; i++)
            {
                rowEnd = i % Math.Sqrt(numberofGenes);
                if(rowEnd == 0)
                {
                    str += " ";
                }
                str += sol[i].ToString() + " ";
            }
            str += " ";
            return str;
        }


        public void Reset()
        {
            IterationID = 0;
            int ThreeTimeSize = populationSize * 3;
            indices = new int[ThreeTimeSize];
            objectiveValue = new double[ThreeTimeSize];
            fitnessValue = new double[ThreeTimeSize];
            Chromosomes = new T[ThreeTimeSize][];
            selectedObjectiveValue = new double[populationSize];
            selectedChromosome = new T[populationSize][];
            for (int r = 0; r < populationSize; r++)
            {
                selectedChromosome[r] = new T[numberofGenes];
            }
            for (int r = 0; r < ThreeTimeSize; r++)
            {
                Chromosomes[r] = new T[numberofGenes];
            }
            mutatedGenes = new bool[populationSize][];
            for (int i = 0; i < populationSize; i++)
            {
                mutatedGenes[i] = new bool[numberofGenes];
            }

            InitializePopulation();

            NumberofCrossedChildren = 0;
            NumberofMutatedChildren = 0;

            if(OptimizationType == GAOptimizationType.Maximization)
            {
                SoFartheBestObjectSolutionValue = double.MinValue;
            }    
            else
            {
                SoFartheBestObjectSolutionValue = double.MaxValue;
            }
            soFartheBestSoltion = new T[numberofGenes];
            UpdateSoFartheBestandIterationtheBest();
            PerformSelectionOperation();
            initialBestObjective = SoFartheBestObjectSolutionValue;
            seriesAverage.Points.Clear();
            seriesIterationtheBest.Points.Clear();
            seriesSoFartheBest.Points.Clear();
        }
        public void RunOneIteration()
        {
            if (IterationID >= IterationLimit) return;
            if(IterationID != 0)
            {
                CopyAllSelectionChromosomes();
            }
            PerformCrossoverOperation();
            PerformMutationOperation();
            UpdateSoFartheBestandIterationtheBest();
            PerformSelectionOperation();
            SeriesUpdate();
            IterationID++;
        }
        public void RunToEnd(int iteration)
        {
            for (int i = 0; i < iteration; i++)
            {
                RunOneIteration();
            }

        }
        void SeriesUpdate()
        {
            SeriesAverage.Points.AddXY(IterationID, iterationAverageObjective);
            SeriesSoFartheBest.Points.AddXY(IterationID, SoFartheBestObjectSolutionValue);
            SeriesIterationtheBest.Points.AddXY(IterationID, iterationBestObjective);
        }
        public enum Optimization
        {
            Minimization = 0,
            Maximization = 1,
            GoalMatching = 2
        }










    }
}
