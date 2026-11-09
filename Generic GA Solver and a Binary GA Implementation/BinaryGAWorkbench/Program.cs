using System;
using System.Globalization;
using System.IO;
using System.Runtime.InteropServices;
using System.Text;
using System.Windows.Forms;
using GALibrary;

namespace BinaryGAWorkbench
{
    internal static class Program
    {
        /// <summary>
        /// Application entry point.
        ///
        /// No arguments        : normal GUI.
        /// --demo              : load the sample matrix, reset and run N generations, then show the GUI.
        /// --demo --headless   : same computation, print the report and exit (no GUI).
        /// --selftest          : run the operator/engine validation suite and exit.
        /// </summary>
        [STAThread]
        private static int Main(string[] args)
        {
            try
            {
                return RunCore(args);
            }
            catch (Exception ex)
            {
                EnsureConsole();
                string text = "Unhandled " + ex.GetType().FullName + ": " + ex.Message
                    + Environment.NewLine + ex.StackTrace;
                Console.WriteLine(text);
                try
                {
                    File.WriteAllText(
                        Path.Combine(AppDomain.CurrentDomain.BaseDirectory, "crash-report.txt"),
                        text,
                        new UTF8Encoding(false));
                }
                catch
                {
                    // Nothing else can be done if the file system is unavailable.
                }
                return 1;
            }
        }

        private static int RunCore(string[] args)
        {
            Options options;
            try
            {
                options = Options.Parse(args, defaultUseBinary: true);
            }
            catch (ArgumentException ex)
            {
                EnsureConsole();
                Console.WriteLine("Argument error: " + ex.Message);
                PrintUsage();
                return 2;
            }

            if (options.Help)
            {
                EnsureConsole();
                PrintUsage();
                return 0;
            }

            if (options.SelfTest)
            {
                EnsureConsole();
                return RunSelfTest(options);
            }

            if (options.Demo && options.Headless)
            {
                EnsureConsole();
                return RunHeadlessDemo(options);
            }

            Application.EnableVisualStyles();
            Application.SetCompatibleTextRenderingDefault(false);

            MainForm form = new MainForm();
            if (options.Demo)
            {
                try
                {
                    CostMatrixProblem sample = LoadSample(options);
                    if (!form.PrepareDemo(sample, options.UseBinary, options.Iterations))
                    {
                        return 1;
                    }
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Could not prepare the demo:" + Environment.NewLine + Environment.NewLine + ex.Message,
                        "Demo error", MessageBoxButtons.OK, MessageBoxIcon.Error);
                    return 1;
                }
            }

            Application.Run(form);
            return 0;
        }

        // -----------------------------------------------------------------
        // Headless modes
        // -----------------------------------------------------------------

        private static int RunHeadlessDemo(Options options)
        {
            try
            {
                CostMatrixProblem problem = LoadSample(options);
                string report = DemoSuite.RunCostMatrixDemos(problem, options.Iterations, options.Seed, null);
                Console.WriteLine(report);
                WriteReport(options.ReportPath, report, "demo-report.txt");
                return 0;
            }
            catch (Exception ex)
            {
                Console.WriteLine("DEMO FAILED: " + ex.GetType().Name + ": " + ex.Message);
                return 1;
            }
        }

        private static int RunSelfTest(Options options)
        {
            SelfTestReport report = GeneticAlgorithmSelfTest.RunAll();
            string text = report.ToText();
            Console.WriteLine(text);
            WriteReport(options.ReportPath, text, "selftest-report.txt");
            return report.Passed ? 0 : 1;
        }

        private static CostMatrixProblem LoadSample(Options options)
        {
            string path = SampleLocator.Find(options.SamplePath, "matrix5.aop");
            if (path == null)
            {
                throw new FileNotFoundException(
                    "No sample matrix was found. Expected samples\\matrix5.aop next to the executable "
                    + "or in a parent directory. Use --sample <path> to select a file.");
            }
            return CostMatrixProblem.Load(path);
        }

        // -----------------------------------------------------------------
        // Console support for a WinExe
        // -----------------------------------------------------------------

        private const int AttachParentProcess = -1;

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool AttachConsole(int processId);

        [DllImport("kernel32.dll", SetLastError = true)]
        private static extern bool AllocConsole();

        private static bool consoleReady;

        private static void EnsureConsole()
        {
            if (consoleReady) return;
            consoleReady = true;
            try
            {
                bool attached = AttachConsole(AttachParentProcess);
                if (!attached)
                {
                    AllocConsole();
                }
                Stream stdout = Console.OpenStandardOutput();
                StreamWriter writer = new StreamWriter(stdout, new UTF8Encoding(false));
                writer.AutoFlush = true;
                Console.SetOut(writer);
                Console.SetError(writer);
            }
            catch
            {
                // When no console is available the report file is still written.
            }
        }

        private static void WriteReport(string requestedPath, string text, string defaultFileName)
        {
            string path = requestedPath;
            if (string.IsNullOrEmpty(path))
            {
                path = Path.Combine(AppDomain.CurrentDomain.BaseDirectory, defaultFileName);
            }
            try
            {
                File.WriteAllText(path, text, new UTF8Encoding(false));
                Console.WriteLine();
                Console.WriteLine("Report written to: " + path);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not write the report: " + ex.Message);
            }
        }

        // -----------------------------------------------------------------
        // Command line
        // -----------------------------------------------------------------

        private sealed class Options
        {
            public bool Demo;
            public bool Headless;
            public bool SelfTest;
            public bool Help;
            public bool UseBinary;
            public int Iterations = 60;
            public int Seed = 2026;
            public string SamplePath;
            public string ReportPath;

            public static Options Parse(string[] args, bool defaultUseBinary)
            {
                Options options = new Options();
                options.UseBinary = defaultUseBinary;
                bool encodingExplicit = false;

                for (int i = 0; i < args.Length; i++)
                {
                    string argument = args[i];
                    if (argument == null) continue;
                    string lower = argument.Trim().ToLowerInvariant();
                    switch (lower)
                    {
                        case "--demo":
                            options.Demo = true;
                            break;
                        case "--headless":
                            options.Headless = true;
                            break;
                        case "--selftest":
                        case "--self-test":
                            options.SelfTest = true;
                            break;
                        case "--binary":
                            options.UseBinary = true;
                            encodingExplicit = true;
                            break;
                        case "--permutation":
                            options.UseBinary = false;
                            encodingExplicit = true;
                            break;
                        case "--help":
                        case "-h":
                        case "/?":
                            options.Help = true;
                            break;
                        case "--run":
                            options.Iterations = ReadInt(args, ref i, "--run", 1, 100000);
                            break;
                        case "--seed":
                            options.Seed = ReadInt(args, ref i, "--seed", int.MinValue, int.MaxValue);
                            break;
                        case "--sample":
                            options.SamplePath = ReadValue(args, ref i, "--sample");
                            break;
                        case "--report":
                            options.ReportPath = ReadValue(args, ref i, "--report");
                            break;
                        default:
                            throw new ArgumentException("Unknown option '" + argument + "'. Use --help for usage.");
                    }
                }

                if (options.Headless && !options.Demo)
                {
                    throw new ArgumentException("--headless requires --demo.");
                }
                if (!encodingExplicit)
                {
                    options.UseBinary = defaultUseBinary;
                }
                return options;
            }

            private static string ReadValue(string[] args, ref int index, string option)
            {
                if (index + 1 >= args.Length)
                {
                    throw new ArgumentException(option + " requires a value.");
                }
                index++;
                return args[index];
            }

            private static int ReadInt(string[] args, ref int index, string option, int minimum, int maximum)
            {
                string text = ReadValue(args, ref index, option);
                int value;
                if (!int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out value)
                    || value < minimum || value > maximum)
                {
                    throw new ArgumentException(option + " requires an integer between "
                        + minimum.ToString(CultureInfo.InvariantCulture) + " and "
                        + maximum.ToString(CultureInfo.InvariantCulture) + ".");
                }
                return value;
            }
        }

        private static void PrintUsage()
        {
            Console.WriteLine("BinaryGAWorkbench - binary-encoded genetic algorithm workbench");
            Console.WriteLine();
            Console.WriteLine("Usage:");
            Console.WriteLine("  BinaryGAWorkbench.exe                          Start the GUI.");
            Console.WriteLine("  BinaryGAWorkbench.exe --demo [--run N]         Load samples\\matrix5.aop,");
            Console.WriteLine("                                                   reset and run N generations in the GUI.");
            Console.WriteLine("  BinaryGAWorkbench.exe --demo --headless [--run N]");
            Console.WriteLine("                                                   Same computation without a window.");
            Console.WriteLine("  BinaryGAWorkbench.exe --selftest               Run the operator/engine test suite.");
            Console.WriteLine();
            Console.WriteLine("Options:");
            Console.WriteLine("  --run N            generations for --demo (default 60).");
            Console.WriteLine("  --sample PATH      cost matrix to load (default samples\\matrix5.aop).");
            Console.WriteLine("  --binary|--permutation  encoding used by the GUI --demo (default binary).");
            Console.WriteLine("  --seed N           deterministic random seed (default 2026).");
            Console.WriteLine("  --report PATH      write the headless report to PATH.");
            Console.WriteLine("  --help             show this help.");
        }
    }

    internal static class SampleLocator
    {
        public static string Find(string explicitPath, string fileName)
        {
            if (!string.IsNullOrEmpty(explicitPath))
            {
                string full = Path.GetFullPath(explicitPath);
                if (!File.Exists(full))
                {
                    throw new FileNotFoundException("Sample file not found: " + full, full);
                }
                return full;
            }

            string directory = AppDomain.CurrentDomain.BaseDirectory;
            for (int depth = 0; depth < 8 && !string.IsNullOrEmpty(directory); depth++)
            {
                string candidate = Path.Combine(directory, "samples", fileName);
                if (File.Exists(candidate)) return candidate;

                candidate = Path.Combine(directory, fileName);
                if (File.Exists(candidate)) return candidate;

                DirectoryInfo parent = Directory.GetParent(directory);
                if (parent == null) break;
                directory = parent.FullName;
            }

            string current = Path.Combine(Environment.CurrentDirectory, "samples", fileName);
            if (File.Exists(current)) return current;
            current = Path.Combine(Environment.CurrentDirectory, fileName);
            if (File.Exists(current)) return current;
            return null;
        }
    }
}


