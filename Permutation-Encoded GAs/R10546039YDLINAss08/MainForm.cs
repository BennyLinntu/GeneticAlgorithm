using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using GALibrary;

namespace R10546039YDLINAss08
{
    public partial class MainForm : Form
    {
        // data fields
        int numberofjobs;
        Random r = new Random();
        GenericSolver<byte> GABinarySolver = null;
        GenericSolver<int> GAPermutationSolver = null;
        
        //functions
        public MainForm()
        {
            InitializeComponent();
        }
        // here are random generate 
        private void btnRandomGenerate_Click(object sender, EventArgs e)
        {
            // clear the dgv and set the dgv
            dgvShow.Rows.Clear();
            dgvShow.Columns.Clear();
            dgvShow.TopLeftHeaderCell.Value = "Job / Machine";

            // record the data
            numberofjobs = int.Parse(tbNumberofJob.Text); // here we record the number of jobs
            double min = double.Parse(tbMinimum.Text);
            double max = double.Parse(tbMaximum.Text);

            // add text
            for (int i = 0; i <= numberofjobs; i++)
            {
                dgvShow.Columns.Add($"Machine{i + 1}", $"Machine{i + 1}");
                dgvShow.Rows.Add();
                dgvShow.Rows[i].HeaderCell.Value = $"Job{i + 1}";
            }

            // add random data
            for (int i = 0; i < dgvShow.Columns.Count; i++)
            {
                for (int j = 0; j < dgvShow.Rows.Count; j++)
                {
                    if(i != j)// if columns are not equal to the rows
                    {
                        double value = Math.Round(min + (r.NextDouble() * (max  - min)), 3);
                        dgvShow[i, j].Value = value; // let it be random
                    }
                    else
                    {
                        dgvShow[i, j].Value = 0.0; // cross point be o
                    }
                }
            }

        }

        // reset the data
        private void btnReset_Click(object sender, EventArgs e)
        {
            // data fields
            int populationSize = int.Parse(tbPopulationSize.Text);
            double crossoverRate = double.Parse(tbCrossoverRate.Text);
            double mutateRate = double.Parse(tbMutationRate.Text);
            double penaltyFactor = double.Parse(tbPenaltyFactor.Text);
            double leastFitnessFraction = double.Parse(tbLeastFitnessFraction.Text);
            GASelectionType selectType = GetSelectionType();
            GAOptimizationType optimizationType = GetGAOptimizationType();
            GAMutationType mutationType = GetBinaryMutationType();
            if(rbBinary.Checked)
            {
                GABinarySolver = new BinaryGA(numberofjobs * numberofjobs, GAOptimizationType.Minimization, GetSetupTimeTotalBInary, GetBinaryCrossoverType());
                GABinarySolver.SelectionType = selectType;
                GABinarySolver.PopulationSize = populationSize;
                GABinarySolver.IterationLimit = 100;
                GABinarySolver.LeastFinessFraction = leastFitnessFraction;
                GABinarySolver.MutationType = mutationType;
                GABinarySolver.OptimizationType = optimizationType;
                GABinarySolver.PenaltyFactor = penaltyFactor;
                GABinarySolver.CrossoverRate = crossoverRate;
                GABinarySolver.MutationRate = mutateRate;
                GABinarySolver.Reset();


            }
            else
            {
                GAPermutationSolver = new PermutationGA(numberofjobs, GAOptimizationType.Minimization, GetSetupTimeTotalPermutation,GetPermutationCrossoverType(), GetPermuatationMutationType());
                GAPermutationSolver.PopulationSize = populationSize;
                GAPermutationSolver.IterationLimit = 100;
                GAPermutationSolver.LeastFinessFraction = leastFitnessFraction;
                GAPermutationSolver.MutationType = mutationType;
                GAPermutationSolver.OptimizationType = optimizationType;
                GAPermutationSolver.PenaltyFactor = penaltyFactor;
                GAPermutationSolver.CrossoverRate = crossoverRate;
                GAPermutationSolver.MutationRate = mutateRate;
                GAPermutationSolver.Reset();
            }
            ClearandAddSeries();
            UpdateNewSolution();

        }

        // get all solution
      

        private void btnGetAllSolution_Click(object sender, EventArgs e)
        {
            int iteration = 100;
            if(rbBinary.Checked)
            {
                GABinarySolver.RunToEnd(iteration);
            }
            else
            {
                GAPermutationSolver.RunToEnd(iteration);
                UpdateNewSolution();
            }
            chtShow.ChartAreas[0].RecalculateAxesScale();
        }

        // Open data
        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (dlgOpen.ShowDialog() != DialogResult.OK) return;
            StreamReader sr = new StreamReader(dlgOpen.FileName);

            double[,] data = null;

            string str = sr.ReadLine();
            numberofjobs = int.Parse(str);

            data = new double[numberofjobs, numberofjobs];

            for (int i = 0; i < numberofjobs; i++)
            {
                str = sr.ReadLine();
                string[] items = str.Split(' ');
                for (int j = 0; j < numberofjobs; j++)
                {
                    data[i, j] = double.Parse(items[j]);
                }
            }
            Read_Problem_at_DGV(data);
            sr.Close();
        }

        // Save Data
        private void saveToolStripMenuItem_Click(object sender, EventArgs e)
        {
            if (dlgSave.ShowDialog() != DialogResult.OK) return;
            numberofjobs = dgvShow.Rows.Count;

            dlgOpen.FileName += ".aop";
            StreamWriter sw = new StreamWriter(dlgSave.FileName);
            sw.WriteLine(numberofjobs);
            for (int i = 0; i < numberofjobs; i++)
            {
                string str = "";
                for (int j = 0; j < numberofjobs; j++)
                {
                    str += dgvShow[j, i].Value.ToString();
                    str += " ";
                }

                sw.WriteLine(str);
            }
            sw.Close();
        }
        //functions
        public double GetSetupTimeTotalBInary(byte[] solution)
        {
            double objective_Value = 0;
            int[] violations = new int[numberofjobs * 2];
            int current_Job_ID = 0; // row count

            for (int row = 0; row < numberofjobs; row++)
            {
                for (int column = 0; column < numberofjobs; column++)
                {
                    objective_Value += Convert.ToDouble(solution[current_Job_ID]) * Convert.ToDouble(dgvShow[row, column].Value);
                    current_Job_ID++;
                }
            }
            return objective_Value;
        }
        public double GetSetupTimeTotalPermutation(int[] solution)
        {
            double total_Time = 0;
            for (int i = 0; i < solution.Length; i++)
            {
                total_Time += Convert.ToDouble(dgvShow[i, solution[i]].Value);
            }
            return total_Time;
        }

        public int[] Return_Chromosomes_Violations(byte[] chrom)
        {
            int[] violations = new int[2 * numberofjobs];
            int counts;

            // row wise
            for (int r = 0; r < numberofjobs; r++)
            {
                counts = 0;
                for (int c = 0; c < numberofjobs; c++)
                    counts += chrom[r * numberofjobs + c];

                violations[r] = Math.Abs(counts - 1);
            }
            // col wise
            for (int c = 0; c < numberofjobs; c++)
            {
                counts = 0;
                for (int r = 0; r < numberofjobs; r++)
                    counts += chrom[r * numberofjobs + c];

                violations[numberofjobs + c] = Math.Abs(counts - 1);
            }
            return violations;
        }
        public void Read_Problem_at_DGV(double[,] data)
        {
            // removed existed data
            dgvShow.Rows.Clear();
            dgvShow.Columns.Clear();

            // adding headers to DGV
            dgvShow.TopLeftHeaderCell.Value = "Job/Machine";
            for (int i = 1; i < numberofjobs + 1; i++)
            {
                dgvShow.Columns.Add($"M{i}", $"M{i}");
                dgvShow.Rows.Add();
                dgvShow.Rows[i - 1].HeaderCell.Value = $"J{i}";
            }

            // adding rnd data to DGV
            for (int col = 0; col < dgvShow.Columns.Count; col++)
            {
                for (int row = 0; row < dgvShow.Rows.Count; row++)
                {
                    dgvShow[col, row].Value = data[row, col];
                }
            }
            dgvShow.RowHeadersWidthSizeMode = DataGridViewRowHeadersWidthSizeMode.AutoSizeToAllHeaders;
        }
        public void ClearandAddSeries()
        {
            chtShow.Series.Clear();
            if(rbBinary.Checked)
            {
                chtShow.Series.Add(GABinarySolver.SeriesAverage);
                chtShow.Series.Add(GABinarySolver.SeriesIterationtheBest);
                chtShow.Series.Add(GABinarySolver.SeriesSoFartheBest);
            }
            else
            {
                chtShow.Series.Add(GAPermutationSolver.SeriesAverage);
                chtShow.Series.Add(GAPermutationSolver.SeriesIterationtheBest);
                chtShow.Series.Add(GAPermutationSolver.SeriesSoFartheBest);
            
            }
        }
        public void UpdateNewSolution()
        {
            int populationSize = int.Parse(tbPopulationSize.Text);
            rtbData.Clear();
            if(rbBinary.Checked)
            {
                rtbData.AppendText($"Best Objective Value:{GABinarySolver.SoFartheBestObjectSolutionValue.ToString()}\n");
                rtbData.AppendText($"Hard Constrain Vilolation\n");
                int[] vio = Return_Chromosomes_Violations(GABinarySolver.SoFartheBestSoltion);
                for (int i = 0; i < numberofjobs; i++)
                {
                    rtbData.AppendText($"Row{i + 1}: {vio[i]}\n");
                }
                for (int i = numberofjobs; i < numberofjobs * 2; i++)
                {
                    rtbData.AppendText($"Col{i + 1 - numberofjobs}: {vio[i]}\n");
                }
                rtbData.AppendText($"The Population\n");
                rtbData.AppendText("Current Iteration: " + GABinarySolver.IterationID.ToString() + "\n");
                rtbData.AppendText("Parents: " + "\n");
                for (int i = 0; i < populationSize ; i++)
                {
                    rtbData.AppendText($"P_{i}: " + GABinarySolver.FlattenSolutiontoString(GABinarySolver.Chromosomes[i]) + "obj: " + GABinarySolver.ObjectiveValue[i] + "\n");
                }
                rtbData.AppendText("Crossovered: " + "\n");
                for (int i = Convert.ToInt32(populationSize); i < populationSize + GABinarySolver.NumberofCrossedChildren; i++)
                {
                    rtbData.AppendText($"C_{i - populationSize}: " + GABinarySolver.FlattenSolutiontoString(GABinarySolver.Chromosomes[i]) + "obj: " + GABinarySolver.ObjectiveValue[i] + "\n");
                }
                rtbData.AppendText("Mutated: " + "\n");
                for (int i = Convert.ToInt32(populationSize + GABinarySolver.NumberofCrossedChildren); i < populationSize + GABinarySolver.NumberofCrossedChildren + GABinarySolver.NumberofMutatedChildren; i++)
                {
                    rtbData.AppendText($"M_{i - Convert.ToInt32(populationSize + GABinarySolver.NumberofCrossedChildren)}: " +GABinarySolver.FlattenSolutiontoString(GABinarySolver.Chromosomes[i]) + "obj: " + GABinarySolver.ObjectiveValue[i] + "\n");
                }
                //rtbData.AppendText($"Shortest Time: "+ GetSetupTimeTotalBInary(GABinarySolver.SoFartheBestObjectSolutionValue).ToString());
            }
            else
            {
               
            }




        }

        public GASelectionType GetSelectionType()
        {
            GASelectionType selectType;
            if(rbStochastic.Checked)
            {
                selectType = GASelectionType.Stochastic;
            }
            else
            {
                selectType = GASelectionType.Deterministic;
            }
            return selectType;
        }
        public GAOptimizationType GetGAOptimizationType()
        {
            if(rbMaximum.Checked)
            {
                return GAOptimizationType.Maximization;
            }
            else
            {
                return GAOptimizationType.Minimization;
            }
        }
        public BinaryCrossoverType GetBinaryCrossoverType()
        {
            if (cbBinaryCrossover.SelectedIndex == 0)
            {
                return BinaryCrossoverType.OnePointCut;
            }
            else if (cbBinaryCrossover.SelectedIndex == 1)
            {
                return BinaryCrossoverType.TwoPointCut;
            }
            else
            {
                return BinaryCrossoverType.NPointCut;
            }
        }
        public GAMutationType GetBinaryMutationType()
        {
            if(cbBinaryMutation.SelectedIndex == 0)
            {
                return GAMutationType.GeneNumberBased;
            }
            else
            {
                return GAMutationType.ChromosomesNumberBased;
            }
        }
        public PermutationCrossoverType GetPermutationCrossoverType()
        {
            if(cbPermutationCrossover.SelectedIndex == 0)
            {
                return PermutationCrossoverType.ParialMappedX;
            }
            else if(cbPermutationCrossover.SelectedIndex == 1)
            {
                return PermutationCrossoverType.OrderX;
            }
            else if(cbPermutationCrossover.SelectedIndex == 2)
            {
                return PermutationCrossoverType.PositionBasedX;
            }
            else if(cbPermutationCrossover.SelectedIndex == 3)
            {
                return PermutationCrossoverType.OrderBasedX;
            }
            else if(cbPermutationCrossover.SelectedIndex == 4)
            {
                return PermutationCrossoverType.CycleX;
            }
            else
            {
                return PermutationCrossoverType.SubtourExchange;
            }
        }
        public PermuatationMutationType GetPermuatationMutationType()
        {
            if (cbPermutationMutation.SelectedIndex == 0)
            {
                return PermuatationMutationType.Inversion;
            }
            else if (cbPermutationMutation.SelectedIndex == 1)
            {
                return PermuatationMutationType.Insertion;
            }
            else if (cbPermutationMutation.SelectedIndex == 2)
            {
                return PermuatationMutationType.Displacement;
            }
            else if (cbPermutationMutation.SelectedIndex == 3)
            {
                return PermuatationMutationType.ReciprocalExchange;
            }
            else
            {
                return PermuatationMutationType.Heuristic;
            }

        }

        private void btnRunOne_Click(object sender, EventArgs e)
        {
            rtbData.Clear();
            if(rbBinary.Checked)
            {
                GABinarySolver.RunOneIteration(); 
            }
            else
            {
                GAPermutationSolver.RunOneIteration();
            }
            UpdateNewSolution();
            chtShow.ChartAreas[0].RecalculateAxesScale();

        }
    }
}
