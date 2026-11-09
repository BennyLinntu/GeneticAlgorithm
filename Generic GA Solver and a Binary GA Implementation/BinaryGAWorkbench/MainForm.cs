using System;
using System.Drawing;
using System.Globalization;
using System.Text;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using GALibrary;

namespace BinaryGAWorkbench
{
    /// <summary>
    /// Binary-encoded GA workbench (selection 07).
    ///
    /// The form owns an <see cref="CostMatrixProblem"/> (the setup-cost matrix)
    /// and creates either a <see cref="BinaryGA"/> or a <see cref="PermutationGA"/>
    /// when the user presses Reset. Every handler validates its input and reports
    /// failures through a message box instead of throwing a null-reference or
    /// format exception.
    /// </summary>
    public partial class MainForm : Form
    {
        private const int DefaultIterationLimit = 100;

        private int numberofjobs;
        private int iterationLimit = DefaultIterationLimit;
        private readonly Random rng = new Random();
        private CostMatrixProblem problem;
        private GenericSolver<byte> binarySolver;
        private GenericSolver<int> permutationSolver;

        public MainForm()
        {
            InitializeComponent();
            Text = "Binary GA Workbench";
            EnsureComboSelection();
            AttachOptionHandlers();
            ConfigureChart();
            UpdateOperatorAvailability();
            ShowWelcome();
        }

        // -----------------------------------------------------------------
        // Startup helpers
        // -----------------------------------------------------------------

        private void ConfigureChart()
        {
            if (chtShow.ChartAreas.Count == 0) return;
            chtShow.ChartAreas[0].AxisX.Title = "Iteration";
            chtShow.ChartAreas[0].AxisY.Title = "Objective value";
            chtShow.ChartAreas[0].AxisX.MajorGrid.LineColor = Color.Gainsboro;
            chtShow.ChartAreas[0].AxisY.MajorGrid.LineColor = Color.Gainsboro;
        }

        private void EnsureComboSelection()
        {
            if (cbBinaryCrossover.SelectedIndex < 0) cbBinaryCrossover.SelectedIndex = 0;
            if (cbBinaryMutation.SelectedIndex < 0) cbBinaryMutation.SelectedIndex = 0;
            if (cbPermutationCrossover.SelectedIndex < 0) cbPermutationCrossover.SelectedIndex = 0;
            if (cbPermutationMutation.SelectedIndex < 0) cbPermutationMutation.SelectedIndex = 0;
        }

        private void AttachOptionHandlers()
        {
            rbBinary.CheckedChanged += OptionChanged;
            rbPermutation.CheckedChanged += OptionChanged;
        }

        private void OptionChanged(object sender, EventArgs e)
        {
            UpdateOperatorAvailability();
        }

        private void UpdateOperatorAvailability()
        {
            bool binary = rbBinary.Checked;
            cbBinaryCrossover.Enabled = binary;
            cbBinaryMutation.Enabled = binary;
            cbPermutationCrossover.Enabled = !binary;
            cbPermutationMutation.Enabled = !binary;
        }

        private void ShowWelcome()
        {
            rtbData.Text =
                "Genetic Algorithm Workbench" + Environment.NewLine +
                Environment.NewLine +
                "1. Press Generate (or File > Open) to create an cost matrix." + Environment.NewLine +
                "2. Choose the binary or permutation encoding and its operators." + Environment.NewLine +
                "3. Press Reset to build generation 0, then Run One Iteration or Get All Solutions." + Environment.NewLine +
                Environment.NewLine +
                "Headless demo: BinaryGAWorkbench.exe --demo --headless --run 60" + Environment.NewLine +
                "Self-test   : BinaryGAWorkbench.exe --selftest" + Environment.NewLine;
        }

        // -----------------------------------------------------------------
        // Public demo API (used by Program --demo)
        // -----------------------------------------------------------------

        public void LoadProblem(CostMatrixProblem value)
        {
            if (value == null) throw new ArgumentNullException(nameof(value));
            problem = value;
            numberofjobs = value.Size;
            FillGrid(value);
            InvalidateSolvers();
        }

        public bool PrepareDemo(CostMatrixProblem value, bool useBinary, int iterations)
        {
            if (value == null) throw new ArgumentNullException(nameof(value));
            iterationLimit = Math.Max(1, iterations);

            LoadProblem(value);

            tbPopulationSize.Text = "30";
            tbCrossoverRate.Text = "0.9";
            tbMutationRate.Text = "0.1";
            tbLeastFitnessFraction.Text = "0.1";
            tbPenaltyFactor.Text = "100";
            rbMinimum.Checked = true;
            if (useBinary)
            {
                rbBinary.Checked = true;
            }
            else
            {
                rbPermutation.Checked = true;
            }

            if (!TryBuildSolver(true))
            {
                return false;
            }

            RunAllIterations(iterations);
            UpdateNewSolution();
            RefreshChart();
            if (tabControl1 != null)
                tabControl1.SelectedTab = tabPage1;
            return true;
        }

        // -----------------------------------------------------------------
        // Grid <-> problem
        // -----------------------------------------------------------------

        private void FillGrid(CostMatrixProblem value)
        {
            dgvShow.SuspendLayout();
            try
            {
                dgvShow.Rows.Clear();
                dgvShow.Columns.Clear();
                dgvShow.TopLeftHeaderCell.Value = "Job / Machine";

                for (int column = 0; column < value.Size; column++)
                {
                    int index = dgvShow.Columns.Add("M" + (column + 1).ToString(CultureInfo.InvariantCulture),
                        "M" + (column + 1).ToString(CultureInfo.InvariantCulture));
                    dgvShow.Columns[index].SortMode = DataGridViewColumnSortMode.NotSortable;
                    dgvShow.Columns[index].Width = 70;
                }

                for (int row = 0; row < value.Size; row++)
                {
                    dgvShow.Rows.Add();
                    dgvShow.Rows[row].HeaderCell.Value = "J" + (row + 1).ToString(CultureInfo.InvariantCulture);
                }

                for (int row = 0; row < value.Size; row++)
                {
                    for (int column = 0; column < value.Size; column++)
                    {
                        dgvShow[column, row].Value = value[row, column].ToString("0.###", CultureInfo.InvariantCulture);
                    }
                }
                dgvShow.RowHeadersWidthSizeMode = DataGridViewRowHeadersWidthSizeMode.AutoSizeToAllHeaders;
            }
            finally
            {
                dgvShow.ResumeLayout();
            }
        }

        private bool TryReadProblemFromGrid(out CostMatrixProblem result, out string error)
        {
            result = null;
            error = null;

            int rows = dgvShow.Rows.Count;
            int columns = dgvShow.Columns.Count;
            if (rows < 2 || columns < 2)
            {
                error = "Generate or open an cost matrix first (at least 2 x 2).";
                return false;
            }
            if (rows != columns)
            {
                error = "The cost matrix must be square: " + rows.ToString(CultureInfo.InvariantCulture)
                    + " rows and " + columns.ToString(CultureInfo.InvariantCulture) + " columns were found.";
                return false;
            }

            double[,] data = new double[rows, columns];
            for (int row = 0; row < rows; row++)
            {
                for (int column = 0; column < columns; column++)
                {
                    object raw = dgvShow[column, row].Value;
                    if (raw == null || raw.ToString().Trim().Length == 0)
                    {
                        error = "Cell (" + (row + 1).ToString(CultureInfo.InvariantCulture) + ", "
                            + (column + 1).ToString(CultureInfo.InvariantCulture) + ") is empty.";
                        return false;
                    }

                    string text = raw.ToString().Trim();
                    double value;
                    bool parsed = double.TryParse(text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out value)
                        || double.TryParse(text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.CurrentCulture, out value);
                    if (!parsed || double.IsNaN(value) || double.IsInfinity(value))
                    {
                        error = "Cell (" + (row + 1).ToString(CultureInfo.InvariantCulture) + ", "
                            + (column + 1).ToString(CultureInfo.InvariantCulture) + ") is not a finite number: '" + text + "'.";
                        return false;
                    }
                    data[row, column] = value;
                }
            }

            try
            {
                result = new CostMatrixProblem(data);
            }
            catch (ArgumentException ex)
            {
                error = ex.Message;
                return false;
            }
            return true;
        }

        private void InvalidateSolvers()
        {
            binarySolver = null;
            permutationSolver = null;
            chtShow.Series.Clear();
        }

        // -----------------------------------------------------------------
        // Parameter validation
        // -----------------------------------------------------------------

        private bool TryReadInt(TextBox box, int minimum, int maximum, string name, out int value)
        {
            value = 0;
            if (!int.TryParse(box.Text.Trim(), NumberStyles.Integer, CultureInfo.InvariantCulture, out value)
                || value < minimum || value > maximum)
            {
                MessageBox.Show(this,
                    name + " must be an integer between " + minimum.ToString(CultureInfo.InvariantCulture)
                    + " and " + maximum.ToString(CultureInfo.InvariantCulture) + ".",
                    "Invalid parameter", MessageBoxButtons.OK, MessageBoxIcon.Warning);
                box.Focus();
                box.SelectAll();
                return false;
            }
            return true;
        }

        private bool TryReadDouble(TextBox box, double minimum, double maximum, string name, out double value)
        {
            value = 0.0;
            string text = box.Text.Trim();
            bool parsed = double.TryParse(text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out value)
                || double.TryParse(text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.CurrentCulture, out value);
            if (!parsed || double.IsNaN(value) || double.IsInfinity(value) || value < minimum || value > maximum)
            {
                MessageBox.Show(this,
                    name + " must be a number between "
                    + minimum.ToString("0.######", CultureInfo.InvariantCulture) + " and "
                    + maximum.ToString("0.######", CultureInfo.InvariantCulture) + ".",
                    "Invalid parameter", MessageBoxButtons.OK, MessageBoxIcon.Warning);
                box.Focus();
                box.SelectAll();
                return false;
            }
            return true;
        }

        private bool TryReadParameters(out int populationSize, out double crossoverRate, out double mutationRate,
            out double penaltyFactor, out double leastFitnessFraction)
        {
            populationSize = 0;
            crossoverRate = 0.0;
            mutationRate = 0.0;
            penaltyFactor = 0.0;
            leastFitnessFraction = 0.0;

            if (!TryReadInt(tbPopulationSize, 2, 5000, "Population size", out populationSize)) return false;
            if (!TryReadDouble(tbCrossoverRate, 0.0, 1.0, "Crossover rate", out crossoverRate)) return false;
            if (!TryReadDouble(tbMutationRate, 0.0, 1.0, "Mutation rate", out mutationRate)) return false;
            if (!TryReadDouble(tbPenaltyFactor, 0.0, 1e12, "Penalty factor", out penaltyFactor)) return false;
            if (!TryReadDouble(tbLeastFitnessFraction, 1e-9, 1.0, "Least fitness fraction", out leastFitnessFraction)) return false;
            return true;
        }

        // -----------------------------------------------------------------
        // Solver construction
        // -----------------------------------------------------------------

        private bool TryBuildSolver(bool showErrors)
        {
            int populationSize;
            double crossoverRate;
            double mutationRate;
            double penaltyFactor;
            double leastFitnessFraction;
            if (!TryReadParameters(out populationSize, out crossoverRate, out mutationRate, out penaltyFactor, out leastFitnessFraction))
            {
                return false;
            }

            CostMatrixProblem gridProblem;
            string error;
            if (!TryReadProblemFromGrid(out gridProblem, out error))
            {
                if (showErrors)
                {
                    MessageBox.Show(this, error, "Cost matrix", MessageBoxButtons.OK, MessageBoxIcon.Warning);
                }
                return false;
            }

            problem = gridProblem;
            numberofjobs = problem.Size;

            try
            {
                if (rbBinary.Checked)
                {
                    BinaryGA solver = new BinaryGA(
                        numberofjobs * numberofjobs,
                        GetGAOptimizationType(),
                        problem.BinaryObjective,
                        GetBinaryCrossoverType());
                    solver.MutationType = GetBinaryMutationType();
                    ApplyCommonSettings(solver, populationSize, crossoverRate, mutationRate, penaltyFactor, leastFitnessFraction);
                    solver.Reset();
                    binarySolver = solver;
                    permutationSolver = null;
                }
                else
                {
                    PermutationGA solver = new PermutationGA(
                        numberofjobs,
                        GetGAOptimizationType(),
                        problem.PermutationObjective,
                        GetPermutationCrossoverType(),
                        GetPermuatationMutationType());
                    ApplyCommonSettings(solver, populationSize, crossoverRate, mutationRate, penaltyFactor, leastFitnessFraction);
                    solver.Reset();
                    permutationSolver = solver;
                    binarySolver = null;
                }
            }
            catch (Exception ex)
            {
                binarySolver = null;
                permutationSolver = null;
                if (showErrors)
                {
                    ShowError("Could not initialise the genetic algorithm.", ex);
                }
                return false;
            }

            ClearandAddSeries();
            UpdateNewSolution();
            return true;
        }

        private void ApplyCommonSettings<T>(
            GenericSolver<T> solver,
            int populationSize,
            double crossoverRate,
            double mutationRate,
            double penaltyFactor,
            double leastFitnessFraction)
        {
            solver.PopulationSize = populationSize;
            solver.CrossoverRate = crossoverRate;
            solver.MutationRate = mutationRate;
            solver.PenaltyFactor = penaltyFactor;
            solver.LeastFinessFraction = leastFitnessFraction;
            solver.IterationLimit = iterationLimit;
            solver.SelectionType = GetSelectionType();
        }

        private void ClearandAddSeries()
        {
            if (chtShow.ChartAreas.Count == 0) return;
            chtShow.Series.Clear();
            if (binarySolver != null)
            {
                AddSeries(binarySolver.SeriesAverage);
                AddSeries(binarySolver.SeriesIterationtheBest);
                AddSeries(binarySolver.SeriesSoFartheBest);
            }
            else if (permutationSolver != null)
            {
                AddSeries(permutationSolver.SeriesAverage);
                AddSeries(permutationSolver.SeriesIterationtheBest);
                AddSeries(permutationSolver.SeriesSoFartheBest);
            }
            chtShow.ChartAreas[0].RecalculateAxesScale();
        }

        private void AddSeries(Series series)
        {
            if (chtShow.ChartAreas.Count > 0)
            {
                series.ChartArea = chtShow.ChartAreas[0].Name;
            }
            chtShow.Series.Add(series);
        }

        // -----------------------------------------------------------------
        // Operators from the UI
        // -----------------------------------------------------------------

        public GASelectionType GetSelectionType()
        {
            return rbStochastic.Checked ? GASelectionType.Stochastic : GASelectionType.Deterministic;
        }

        public GAOptimizationType GetGAOptimizationType()
        {
            return rbMaximum.Checked ? GAOptimizationType.Maximization : GAOptimizationType.Minimization;
        }

        public BinaryCrossoverType GetBinaryCrossoverType()
        {
            switch (Math.Max(0, cbBinaryCrossover.SelectedIndex))
            {
                case 0: return BinaryCrossoverType.OnePointCut;
                case 1: return BinaryCrossoverType.TwoPointCut;
                default: return BinaryCrossoverType.NPointCut;
            }
        }

        public GAMutationType GetBinaryMutationType()
        {
            return Math.Max(0, cbBinaryMutation.SelectedIndex) == 0
                ? GAMutationType.GeneNumberBased
                : GAMutationType.ChromosomesNumberBased;
        }

        public PermutationCrossoverType GetPermutationCrossoverType()
        {
            switch (Math.Max(0, cbPermutationCrossover.SelectedIndex))
            {
                case 0: return PermutationCrossoverType.ParialMappedX;
                case 1: return PermutationCrossoverType.OrderX;
                case 2: return PermutationCrossoverType.PositionBasedX;
                case 3: return PermutationCrossoverType.OrderBasedX;
                case 4: return PermutationCrossoverType.CycleX;
                default: return PermutationCrossoverType.SubtourExchange;
            }
        }

        public PermuatationMutationType GetPermuatationMutationType()
        {
            switch (Math.Max(0, cbPermutationMutation.SelectedIndex))
            {
                case 0: return PermuatationMutationType.Inversion;
                case 1: return PermuatationMutationType.Insertion;
                case 2: return PermuatationMutationType.Displacement;
                case 3: return PermuatationMutationType.ReciprocalExchange;
                default: return PermuatationMutationType.Heuristic;
            }
        }

        // -----------------------------------------------------------------
        // Button handlers
        // -----------------------------------------------------------------

        private void btnRandomGenerate_Click(object sender, EventArgs e)
        {
            int size;
            double minimum;
            double maximum;
            if (!TryReadInt(tbNumberofJob, 2, 50, "Number of jobs / machines", out size)) return;
            if (!TryReadDouble(tbMinimum, -1e9, 1e9, "Minimum", out minimum)) return;
            if (!TryReadDouble(tbMaximum, -1e9, 1e9, "Maximum", out maximum)) return;
            if (minimum > maximum)
            {
                MessageBox.Show(this, "Minimum must not be larger than Maximum.", "Invalid range",
                    MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return;
            }

            try
            {
                problem = CostMatrixProblem.Generate(size, minimum, maximum, rng);
                numberofjobs = size;
                FillGrid(problem);
                InvalidateSolvers();
                rtbData.Text = "Generated a " + size.ToString(CultureInfo.InvariantCulture)
                    + " x " + size.ToString(CultureInfo.InvariantCulture)
                    + " matrix in [" + minimum.ToString("0.###", CultureInfo.InvariantCulture)
                    + ", " + maximum.ToString("0.###", CultureInfo.InvariantCulture) + "]."
                    + Environment.NewLine + "Press Reset to build generation 0.";
            }
            catch (Exception ex)
            {
                ShowError("Could not generate the matrix.", ex);
            }
        }

        private void btnReset_Click(object sender, EventArgs e)
        {
            TryBuildSolver(true);
        }

        private void btnRunOne_Click(object sender, EventArgs e)
        {
            if (!HasSolver())
            {
                PromptForReset();
                return;
            }

            RunWithWaitCursor(delegate
            {
                if (binarySolver != null) binarySolver.RunOneIteration();
                else if (permutationSolver != null) permutationSolver.RunOneIteration();
            });
            UpdateNewSolution();
            RefreshChart();
        }

        private void btnGetAllSolution_Click(object sender, EventArgs e)
        {
            if (!HasSolver())
            {
                PromptForReset();
                return;
            }

            RunWithWaitCursor(delegate { RunAllIterations(iterationLimit); });
            UpdateNewSolution();
            RefreshChart();
        }

        private void RunAllIterations(int iterations)
        {
            if (binarySolver != null)
            {
                binarySolver.RunToEnd(iterations);
            }
            else if (permutationSolver != null)
            {
                permutationSolver.RunToEnd(iterations);
            }
        }

        private bool HasSolver()
        {
            return binarySolver != null || permutationSolver != null;
        }

        private void PromptForReset()
        {
            MessageBox.Show(this,
                "No genetic algorithm has been initialised yet. Generate or open a matrix and press Reset first.",
                "Reset required", MessageBoxButtons.OK, MessageBoxIcon.Information);
        }

        private void RunWithWaitCursor(Action action)
        {
            System.Windows.Forms.Cursor previous = this.Cursor;
            this.Cursor = Cursors.WaitCursor;
            try
            {
                action();
            }
            catch (Exception ex)
            {
                ShowError("The genetic algorithm failed.", ex);
            }
            finally
            {
                this.Cursor = previous;
            }
        }

        private void RefreshChart()
        {
            if (chtShow.ChartAreas.Count > 0)
            {
                chtShow.ChartAreas[0].RecalculateAxesScale();
            }
        }

        // -----------------------------------------------------------------
        // File menu
        // -----------------------------------------------------------------

        private void openToolStripMenuItem_Click(object sender, EventArgs e)
        {
            dlgOpen.Filter = "Cost-matrix problem (*.aop)|*.aop|Text files (*.txt)|*.txt|All files (*.*)|*.*";
            dlgOpen.Title = "Open cost matrix";
            if (dlgOpen.ShowDialog(this) != DialogResult.OK) return;

            try
            {
                problem = CostMatrixProblem.Load(dlgOpen.FileName);
                numberofjobs = problem.Size;
                FillGrid(problem);
                InvalidateSolvers();
                rtbData.Text = "Loaded " + dlgOpen.FileName + Environment.NewLine
                    + "Matrix size: " + numberofjobs.ToString(CultureInfo.InvariantCulture)
                    + " x " + numberofjobs.ToString(CultureInfo.InvariantCulture)
                    + Environment.NewLine + "Press Reset to build generation 0.";
            }
            catch (Exception ex)
            {
                ShowError("Could not open the cost-matrix file.", ex);
            }
        }

        private void saveToolStripMenuItem_Click(object sender, EventArgs e)
        {
            CostMatrixProblem current;
            string error;
            if (!TryReadProblemFromGrid(out current, out error))
            {
                MessageBox.Show(this, error, "Cost matrix", MessageBoxButtons.OK, MessageBoxIcon.Warning);
                return;
            }

            if (!HasSolver())
            {
                problem = current;
            }
            numberofjobs = current.Size;

            dlgSave.Filter = "Cost-matrix problem (*.aop)|*.aop|Text files (*.txt)|*.txt|All files (*.*)|*.*";
            dlgSave.DefaultExt = "aop";
            dlgSave.AddExtension = true;
            dlgSave.Title = "Save cost matrix";
            dlgSave.FileName = "selection" + numberofjobs.ToString(CultureInfo.InvariantCulture) + ".aop";
            if (dlgSave.ShowDialog(this) != DialogResult.OK) return;

            try
            {
                current.Save(dlgSave.FileName);
                rtbData.Text = "Saved " + numberofjobs.ToString(CultureInfo.InvariantCulture)
                    + " x " + numberofjobs.ToString(CultureInfo.InvariantCulture)
                    + " matrix to " + dlgSave.FileName + ".";
            }
            catch (Exception ex)
            {
                ShowError("Could not save the cost-matrix file.", ex);
            }
        }

        // -----------------------------------------------------------------
        // Results pane
        // -----------------------------------------------------------------

        private void UpdateNewSolution()
        {
            if (!HasSolver())
            {
                ShowWelcome();
                return;
            }

            StringBuilder builder = new StringBuilder();
            if (binarySolver != null)
            {
                AppendBinaryReport(builder, binarySolver);
            }
            else
            {
                AppendPermutationReport(builder, permutationSolver);
            }
            rtbData.Text = builder.ToString();
        }

        private void AppendBinaryReport(StringBuilder builder, GenericSolver<byte> solver)
        {
            builder.AppendLine("Encoding          : Binary (n x n cost matrix)");
            builder.AppendLine("Optimization      : " + solver.OptimizationType);
            builder.AppendLine("Selection         : " + solver.SelectionType);
            builder.AppendLine("Crossover         : " + GetBinaryCrossoverType() + " (rate " + solver.CrossoverRate.ToString("0.###", CultureInfo.InvariantCulture) + ")");
            builder.AppendLine("Mutation          : " + GetBinaryMutationType() + " (rate " + solver.MutationRate.ToString("0.###", CultureInfo.InvariantCulture) + ")");
            builder.AppendLine("Population / limit: " + solver.PopulationSize.ToString(CultureInfo.InvariantCulture)
                + " / " + solver.IterationLimit.ToString(CultureInfo.InvariantCulture));
            builder.AppendLine("Iteration         : " + solver.IterationID.ToString(CultureInfo.InvariantCulture)
                + " / " + solver.IterationLimit.ToString(CultureInfo.InvariantCulture));
            builder.AppendLine("Best penalised    : " + solver.SoFartheBestObjectSolutionValue.ToString("0.####", CultureInfo.InvariantCulture));
            builder.AppendLine("Iteration best    : " + solver.IterationBestObjective.ToString("0.####", CultureInfo.InvariantCulture));
            builder.AppendLine("Iteration average : " + solver.IterationAverageObjective.ToString("0.####", CultureInfo.InvariantCulture));
            builder.AppendLine("Penalty factor    : " + solver.PenaltyFactor.ToString("0.###", CultureInfo.InvariantCulture));

            byte[] best = solver.SoFartheBestSoltion;
            if (best != null && problem != null)
            {
                double raw = problem.BinaryObjective(best);
                int[] violations = problem.BinaryViolations(best);
                int totalViolations = 0;
                for (int i = 0; i < violations.Length; i++) totalViolations += violations[i];

                builder.AppendLine("Best raw cost     : " + raw.ToString("0.####", CultureInfo.InvariantCulture));
                builder.AppendLine("Violations        : " + totalViolations.ToString(CultureInfo.InvariantCulture)
                    + " (rows + columns, sum of absolute deviations)");

                builder.AppendLine();
                builder.AppendLine("Constraint detail:");
                for (int row = 0; row < numberofjobs; row++)
                {
                    builder.AppendLine("  Row J" + (row + 1).ToString(CultureInfo.InvariantCulture) + ": " + violations[row].ToString(CultureInfo.InvariantCulture));
                }
                for (int column = 0; column < numberofjobs; column++)
                {
                    builder.AppendLine("  Col M" + (column + 1).ToString(CultureInfo.InvariantCulture) + ": "
                        + violations[numberofjobs + column].ToString(CultureInfo.InvariantCulture));
                }

                builder.AppendLine();
                builder.AppendLine("Best selection (1 = selected):");
                for (int row = 0; row < numberofjobs; row++)
                {
                    StringBuilder line = new StringBuilder("  J" + (row + 1).ToString(CultureInfo.InvariantCulture) + ": ");
                    int assigned = 0;
                    for (int column = 0; column < numberofjobs; column++)
                    {
                        byte gene = best[row * numberofjobs + column];
                        line.Append(gene == 0 ? "0 " : "1 ");
                        if (gene != 0) assigned = column + 1;
                    }
                    line.Append(" -> M" + assigned.ToString(CultureInfo.InvariantCulture));
                    builder.AppendLine(line.ToString());
                }
                builder.AppendLine("Best chromosome   : " + solver.FlattenSolutiontoString(best));
            }

            AppendPool(builder, solver);
        }

        private void AppendPermutationReport(StringBuilder builder, GenericSolver<int> solver)
        {
            builder.AppendLine("Encoding          : Permutation");
            builder.AppendLine("Optimization      : " + solver.OptimizationType);
            builder.AppendLine("Selection         : " + solver.SelectionType);
            builder.AppendLine("Crossover         : " + GetPermutationCrossoverType() + " (rate " + solver.CrossoverRate.ToString("0.###", CultureInfo.InvariantCulture) + ")");
            builder.AppendLine("Mutation          : " + GetPermuatationMutationType() + " (rate " + solver.MutationRate.ToString("0.###", CultureInfo.InvariantCulture) + ")");
            builder.AppendLine("Population / limit: " + solver.PopulationSize.ToString(CultureInfo.InvariantCulture)
                + " / " + solver.IterationLimit.ToString(CultureInfo.InvariantCulture));
            builder.AppendLine("Iteration         : " + solver.IterationID.ToString(CultureInfo.InvariantCulture)
                + " / " + solver.IterationLimit.ToString(CultureInfo.InvariantCulture));
            builder.AppendLine("Best setup cost   : " + solver.SoFartheBestObjectSolutionValue.ToString("0.####", CultureInfo.InvariantCulture));
            builder.AppendLine("Iteration best    : " + solver.IterationBestObjective.ToString("0.####", CultureInfo.InvariantCulture));
            builder.AppendLine("Iteration average : " + solver.IterationAverageObjective.ToString("0.####", CultureInfo.InvariantCulture));

            int[] best = solver.SoFartheBestSoltion;
            if (best != null && problem != null)
            {
                builder.AppendLine("Best permutation  : " + FormatPermutation(best));
                builder.AppendLine();
                builder.AppendLine("Selection:");
                try
                {
                    for (int job = 0; job < best.Length; job++)
                    {
                        builder.AppendLine("  Job J" + (job + 1).ToString(CultureInfo.InvariantCulture)
                            + " -> Machine M" + (best[job] + 1).ToString(CultureInfo.InvariantCulture)
                            + "   cost " + problem[job, best[job]].ToString("0.###", CultureInfo.InvariantCulture));
                    }
                }
                catch (ArgumentException ex)
                {
                    builder.AppendLine("  The best chromosome is not a valid permutation: " + ex.Message);
                }
            }

            AppendPool(builder, solver);
        }

        private void AppendPool(StringBuilder builder, GenericSolver<byte> solver)
        {
            builder.AppendLine();
            builder.AppendLine("Active pool rows  : " + solver.ActivePoolSize.ToString(CultureInfo.InvariantCulture)
                + " (N + C + M, capacity " + (solver.PopulationSize * 3).ToString(CultureInfo.InvariantCulture) + ")");
            builder.AppendLine("Parents:");
            int shown = 0;
            for (int i = 0; i < solver.PopulationSize && shown < 30; i++, shown++)
            {
                AppendRow(builder, "  P_" + i.ToString(CultureInfo.InvariantCulture), solver, i);
            }
            builder.AppendLine("Crossover children:");
            for (int i = solver.PopulationSize; i < solver.PopulationSize + solver.NumberofCrossedChildren && shown < 40; i++, shown++)
            {
                AppendRow(builder, "  C_" + (i - solver.PopulationSize).ToString(CultureInfo.InvariantCulture), solver, i);
            }
            builder.AppendLine("Mutation children:");
            int mutationStart = solver.PopulationSize + solver.NumberofCrossedChildren;
            for (int i = mutationStart; i < mutationStart + solver.NumberofMutatedChildren && shown < 50; i++, shown++)
            {
                AppendRow(builder, "  M_" + (i - mutationStart).ToString(CultureInfo.InvariantCulture), solver, i);
            }
            if (shown >= 50)
            {
                builder.AppendLine("  ... display truncated; the chart and best values use the full pool.");
            }
        }

        private void AppendPool(StringBuilder builder, GenericSolver<int> solver)
        {
            builder.AppendLine();
            builder.AppendLine("Active pool rows  : " + solver.ActivePoolSize.ToString(CultureInfo.InvariantCulture)
                + " (N + C + M, capacity " + (solver.PopulationSize * 3).ToString(CultureInfo.InvariantCulture) + ")");
            builder.AppendLine("Parents:");
            int shown = 0;
            for (int i = 0; i < solver.PopulationSize && shown < 30; i++, shown++)
            {
                AppendRow(builder, "  P_" + i.ToString(CultureInfo.InvariantCulture), solver, i);
            }
            builder.AppendLine("Crossover children:");
            for (int i = solver.PopulationSize; i < solver.PopulationSize + solver.NumberofCrossedChildren && shown < 40; i++, shown++)
            {
                AppendRow(builder, "  C_" + (i - solver.PopulationSize).ToString(CultureInfo.InvariantCulture), solver, i);
            }
            builder.AppendLine("Mutation children:");
            int mutationStart = solver.PopulationSize + solver.NumberofCrossedChildren;
            for (int i = mutationStart; i < mutationStart + solver.NumberofMutatedChildren && shown < 50; i++, shown++)
            {
                AppendRow(builder, "  M_" + (i - mutationStart).ToString(CultureInfo.InvariantCulture), solver, i);
            }
            if (shown >= 50)
            {
                builder.AppendLine("  ... display truncated; the chart and best values use the full pool.");
            }
        }

        private static void AppendRow<T>(StringBuilder builder, string label, GenericSolver<T> solver, int index)
        {
            if (index < 0 || index >= solver.ActivePoolSize) return;
            builder.AppendLine(label + ": " + solver.FlattenSolutiontoString(solver.Chromosomes[index])
                + " obj=" + solver.ObjectiveValue[index].ToString("0.####", CultureInfo.InvariantCulture));
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

        private void ShowError(string context, Exception ex)
        {
            MessageBox.Show(this,
                context + Environment.NewLine + Environment.NewLine + ex.Message,
                "Error", MessageBoxButtons.OK, MessageBoxIcon.Error);
        }
    }
}



