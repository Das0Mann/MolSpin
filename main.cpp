//////////////////////////////////////////////////////////////////////////////
// MolSpin - Main program file
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
// Adjust this number to fit your system
#define MAX_THREADS 48
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "Settings.h"
#include "MSDParser.h"
#include "RunSection.h"
#include "FileReader.h"
#include <fstream>
#include <unistd.h>

//////////////////////////////////////////////////////////////////////////////
// #ifdef USE_OPENBLAS
extern "C" void openblas_set_num_threads(int);
// #endif
// #ifdef USE_OPENMP
extern "C" void omp_set_num_threads(int);
// #endif
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	const std::string MolSpin_version = "v2.0";
	const std::string hline = "# ---------------------------------------------------------------------------------";

	// Check for proper input
	if (argc < 2)
	{
		std::cout << "Please specify an input file:" << std::endl;
		std::cout << "  molspin [inputfile]" << std::endl;
		return 0;
	}

	// Create a temporary string to hold the contents of "argv"
	std::string strargv = argv[1];

	// Use -h or --help to get a help message
	if (argc == 2 && (strargv.compare("-h") == 0 || strargv.compare("--help") == 0))
	{
		std::cout << "Usage:\n  molspin [options] [inputfile]\n"
				  << std::endl;
		std::cout << "Options:" << std::endl;
		std::cout << "    -a\n    --append" << std::endl;
		std::cout << "             File output is appended to existing files (instead of overwriting)." << std::endl;
		std::cout << "             Example: molspin -a myfile.msd" << std::endl;
		std::cout << "    -c\n    --checkpoint" << std::endl;
		std::cout << "             Skip all tasks in the runsection until the specified checkpoint is found." << std::endl;
		std::cout << "             This only happens for the first step." << std::endl;
		std::cout << "             Example: molspin -c my_taskname myfile.msd" << std::endl;
		std::cout << "    -d\n    --defines" << std::endl;
		std::cout << "             Show all defined directives and their values, and included files." << std::endl;
		std::cout << "             Example: molspin -d myfile.msd" << std::endl;
		std::cout << "    -h\n    --help" << std::endl;
		std::cout << "             Shows this message if no input file is specified." << std::endl;
		std::cout << "             Example: molspin -h" << std::endl;
		std::cout << "    -n\n    --steps" << std::endl;
		std::cout << "             Specify how many steps should be calculated before you are notified." << std::endl;
		std::cout << "             Example: molspin -n 5 myfile.msd" << std::endl;
		std::cout << "    -o\n    --objects" << std::endl;
		std::cout << "             Show the objects that were read from the input files." << std::endl;
		std::cout << "             Example: molspin -o myfile.msd" << std::endl;
		std::cout << "    -os\n    --objects-states" << std::endl;
		std::cout << "             Same as -o but with full information about State objects." << std::endl;
		std::cout << "             Example: molspin -os myfile.msd" << std::endl;
		std::cout << "    -p\n    --threads" << std::endl;
		std::cout << "             Specify the number of threads/processor cores to use." << std::endl;
		std::cout << "             Example: molspin -p 24 myfile.msd" << std::endl;
		std::cout << "    -r\n    --first-step" << std::endl;
		std::cout << "             Specify which step to start from (if you don't want to start from step 0)." << std::endl;
		std::cout << "             Example: molspin -r 5 myfile.msd" << std::endl;
		std::cout << "    -l\n    --step-limit" << std::endl;
		std::cout << "             Limits the number of steps to run." << std::endl;
		std::cout << "             Example: molspin -l 5 myfile.msd" << std::endl;
		std::cout << "    -s\n    --silent" << std::endl;
		std::cout << "             Minimize the text output." << std::endl;
		std::cout << "             Example: molspin -s myfile.msd" << std::endl;
		std::cout << "    -t\n    --action-targets" << std::endl;
		std::cout << "             Prints all the ActionTargets that can be used." << std::endl;
		std::cout << "             Example: molspin -t myfile.msd" << std::endl;
		std::cout << "    -z\n    --no-calc" << std::endl;
		std::cout << "             Skip calculations." << std::endl;
		std::cout << "             Example: molspin -z myfile.msd" << std::endl;
		std::cout << "\n"
				  << std::endl;
		std::cout << "    Note: You can restart a stopped calculation using:  " << std::endl;
		std::cout << "    molspin -r <step number> -c <task name> -a myfile.msd" << std::endl;
		return 0;
	}

	std::cout << hline << std::endl;
	std::cout << "# 	        &&&&&&&                                                                                                                " << std::endl;
	std::cout << "#      &&&&&&(      *&&&&&&                                                                                                          " << std::endl;
	std::cout << "#      &&&(            ,&&&      &&&&           &&&&                    &&&#      &&&&&&&/                        &&&                " << std::endl;
	std::cout << "#        &&&&&&&&&&&&&&&&        &&&&&&       &&&&&&                    &&&#   /&&&&& &&&&&&                      &&&                " << std::endl;
	std::cout << "#           &&&    &&&#          &&&&&&&    #&&&&&&&       &&&&&.       &&&#   &&&/             &&&&&&&&&&&%             &&  &&&&&   " << std::endl;
	std::cout << "#    ,&&&&&&  &&&&&&  &&&&&&&    &&&  &&&& &&&& &&&&    &&&&&%&&&&&&    &&&#   &&&&&&&          &&&&&&&#&&&&&%    &&&    &&&&&&&&&&&&" << std::endl;
	std::cout << "# &&&.                       &&& &&&    &&&&&   &&&&   &&&        &&&   &&&#      (&&&&&&&&&    &&&&       &&&&   &&&    &&&      &&&" << std::endl;
	std::cout << "#    &&&&&&&&&&&&&&&&&&&&&&&&/   &&&     &&     &&&&  &&&&        &&&   &&&#             %&&&   &&&        &&&&   &&&    &&&      &&&" << std::endl;
	std::cout << "#        &&&&      &&&&          &&&            &&&&   &&&&      &&&&   &&&#   &&&       &&&&   &&&&      &&&&    &&&    &&&      &&&" << std::endl;
	std::cout << "#       &&&            &&&,      &&&            &&&&     &&&&&&&&&&     &&&#    &&&&&&&&&&&&    &&& &&&&&&&&&     &&&    &&&      &&&" << std::endl;
	std::cout << "#      &&/               &&&                                                                    &&&                                  " << std::endl;
	std::cout << "#        &&&&&&&&&&&&&&&&                                                                       &&&      							   " << std::endl;
	std::cout << hline << std::endl;
	std::cout << "# Molecular Spin Dynamics " << MolSpin_version << std::endl;
	std::cout << "# " << std::endl;
	std::cout << "# Developed 2017-2019 by Claus Nielsen* and 2021-2024 by Luca Gerhards**." << std::endl;
	std::cout << "# (c) Quantum Biology and Computational Physics Group," << std::endl;
	std::cout << "# *:  University of Southern Denmark." << std::endl;
	std::cout << "# **: Carl von Ossietzky University of Oldenburg." << std::endl;
	std::cout << "# For more information see www.molspin.eu" << std::endl;
	std::cout << hline << std::endl;

	// Program options to be set on the commandline
	bool silentMode = false;
	bool printActionTargets = false;
	bool noCalculations = false;
	bool appendMode = false;
	bool hasCheckpoint = false;
	bool showDefines = false;
	bool showObjects = false;
	bool showFullObjects = false;
	unsigned int reportSteps = 1;
	unsigned int firstStep = 1;
	unsigned int stepLimit = 0;
	std::string checkpoint = "";

	// -----------------------------------------------------
	// START Parsing of commandline options
	// -----------------------------------------------------
	// If some options are provided

	std::cout << "# Recognized commandline options will be listed in this section." << std::endl;
	if (argc > 2)
	{
		// Loop through all of them
		for (int i = 1; i < argc - 1; i++)
		{
			// Set the temporary string to the current string in "argv"
			strargv = argv[i];

			// Number of threads specified
			if (strargv.compare("-p") == 0 || strargv.compare("--threads") == 0)
			{
				// Check whether a number is specified
				if (argc - 2 <= i)
				{
					std::cout << "# - Warning: Either number-of-threads specification or inputfile specification is missing!" << std::endl;
					std::cout << "#   Example of usage of -p/--threads:\n#   molspin -p 24 example.msd" << std::endl;
					return 0;
				}

				// Set the number of threads
				try
				{
					int threads = std::stoi(argv[i + 1]);

					if (threads >= 1 && threads <= MAX_THREADS)
					{
						std::cout << "# - Number of threads set to " << threads << "." << std::endl;
						// #ifdef USE_OPENBLAS
						openblas_set_num_threads(threads);
						// #endif
						// #ifdef USE_OPENMP
						omp_set_num_threads(threads);
						// #endif
					}
					else
					{
						std::cout << "# - Could not set number of threads to " << threads << "! Please specify a number from 1 to " << MAX_THREADS << "." << std::endl;
					}
				}
				catch (const std::exception &) // Catch any conversion errors
				{
					std::cout << "# - Could not set number of threads to " << argv[i + 1] << "! Please specify a valid number from 1 to " << MAX_THREADS << "." << std::endl;
				}

				// Don't try to parse the number of threads as a commandline parameter
				i++;
			}
			else if (strargv.compare("-s") == 0 || strargv.compare("--silent") == 0)
			{
				silentMode = true;
				std::cout << "# - Silent-mode." << std::endl;
			}
			else if (strargv.compare("-a") == 0 || strargv.compare("--append") == 0)
			{
				appendMode = true;
				std::cout << "# - Appending file output (instead of overwriting)." << std::endl;
			}
			else if (strargv.compare("-d") == 0 || strargv.compare("--defines") == 0)
			{
				// No need to show this message if "showDefines" was already invoked
				if (!showDefines)
				{
					showDefines = true;
					std::cout << "# - Showing defined names and values, and included files." << std::endl;
				}

				// Check whether a defined name was specified
				if (argc - 2 > i)
				{
					// Get the name and value for the definition
					std::string defineName(argv[i + 1]);
					std::string defineValue = "";

					// Check whether there is a value specified
					if (argc - 3 > i)
					{
						defineValue = std::string(argv[i + 2]);
					}

					// Check that the name was not just the next commandline parameter
					if (!defineName.empty() && defineName[0] != '-')
					{
						// Don't try to parse the defineName as a commandline parameter
						i++;

						// Check that the value was not just the next commandline parameter
						if (!defineValue.empty())
						{
							if (defineValue[0] == '-')
								defineValue = "";
							else
								i++; // Don't try to parse the defineValue as a commandline parameter
						}

						// Add the definition
						MSDParser::FileReader::AddDefinition(defineName, defineValue);
					}
				}
			}
			else if (strargv.compare("-o") == 0 || strargv.compare("--objects") == 0)
			{
				showObjects = true;
				std::cout << "# - Printing objects read from the input files." << std::endl;
			}
			else if (strargv.compare("-os") == 0 || strargv.compare("--objects-states") == 0)
			{
				showObjects = true;
				showFullObjects = true;
				std::cout << "# - Printing objects read from the input files with full state information." << std::endl;
			}
			else if (strargv.compare("-t") == 0 || strargv.compare("--action-targets") == 0)
			{
				printActionTargets = true;
				std::cout << "# - Printing ActionTargets." << std::endl;
			}
			else if (strargv.compare("-z") == 0 || strargv.compare("--no-calc") == 0)
			{
				noCalculations = true;
				std::cout << "# - Skipping calculations." << std::endl;
			}
			else if (strargv.compare("-n") == 0 || strargv.compare("--steps") == 0)
			{
				// Check whether a number is specified
				if (argc - 2 <= i)
				{
					std::cout << "# - Warning: Either steps-done notification frequency or inputfile specification is missing!" << std::endl;
					std::cout << "#   Example of usage of -n/--steps:\n#   molspin -n 5 example.msd" << std::endl;
					return 0;
				}

				// Specify how often (after how many steps) "step done" is printed
				try
				{
					reportSteps = std::stoi(argv[i + 1]);
				}
				catch (const std::exception &) // Catch any conversion errors
				{
					reportSteps = 1;
					std::cout << "# - Could not set steps-done notification frequency to \"" << argv[i + 1] << "\"." << std::endl;
					std::cout << "#   Example of usage of -n/--steps:\n#   molspin -n 10 example.msd" << std::endl;
				}
				std::cout << "# - Steps-done notification frequency set to " << reportSteps << "." << std::endl;

				// Don't try to parse the number of steps as a commandline parameter
				i++;
			}
			else if (strargv.compare("-r") == 0 || strargv.compare("--first-step") == 0)
			{
				// Check whether a number is specified
				if (argc - 2 <= i)
				{
					std::cout << "# - Warning: Either first-step or inputfile specification is missing!" << std::endl;
					std::cout << "#   Example of usage of -r/--firststep:\n#   molspin -r 5 example.msd" << std::endl;
					return 0;
				}

				// Specify the which step number should be the first
				try
				{
					firstStep = std::stoi(argv[i + 1]);
					if (!std::isfinite(firstStep))
					{
						firstStep = 0;
						std::cout << "# - Could not set first-step to \"" << argv[i + 1] << "\"." << std::endl;
						std::cout << "#   Example of usage of -r/--first-step:\n#   molspin -r 10 example.msd" << std::endl;
					}
					std::cout << "# - First step set to " << firstStep << "." << std::endl;
				}
				catch (const std::exception &) // Catch any conversion errors
				{
					std::cout << "# - Could not set first-step to \"" << argv[i + 1] << "\"." << std::endl;
					std::cout << "#   Example of usage of -r/--first-step:\n#   molspin -r 10 example.msd" << std::endl;
				}

				// Don't try to parse the restart step number as a commandline parameter
				i++;
			}
			else if (strargv.compare("-l") == 0 || strargv.compare("--step-limit") == 0)
			{
				// Check whether a number is specified
				if (argc - 2 <= i)
				{
					std::cout << "# - Warning: Either step-limit or inputfile specification is missing!" << std::endl;
					std::cout << "#   Example of usage of -l/--step-limit:\n#   molspin -l 5 example.msd" << std::endl;
					return 0;
				}

				// Specify the maximum number of steps
				try
				{
					stepLimit = std::stoi(argv[i + 1]);
					if (!std::isfinite(stepLimit))
					{
						stepLimit = 0;
						std::cout << "# - Could not set step-limit to \"" << argv[i + 1] << "\"." << std::endl;
						std::cout << "#   Example of usage of -l/--step-limit:\n#   molspin -l 10 example.msd" << std::endl;
					}
					std::cout << "# - Number of steps limited to " << stepLimit << "." << std::endl;
				}
				catch (const std::exception &) // Catch any conversion errors
				{
					std::cout << "# - Could not set step-limit to \"" << argv[i + 1] << "\"." << std::endl;
					std::cout << "#   Example of usage of -l/--step-limit:\n#   molspin -l 10 example.msd" << std::endl;
				}

				// Don't try to parse the step limit as a commandline parameter
				i++;
			}
			else if (strargv.compare("-c") == 0 || strargv.compare("--checkpoint") == 0)
			{
				// Check whether a number is specified
				if (argc - 2 <= i)
				{
					std::cout << "# - Warning: Either checkpoint (task) name or inputfile specification is missing!" << std::endl;
					std::cout << "#   Example of usage of -c/--checkpoint:\n#   molspin -c my_taskname example.msd" << std::endl;
					return 0;
				}

				// Specify how often (after how many steps) "step done" is printed
				checkpoint = argv[i + 1];
				hasCheckpoint = true;
				std::cout << "# - Checkpoint set to \"" << checkpoint << "\"." << std::endl;

				// Don't try to parse the checkpoint name as a commandline parameter
				i++;
			}
			else
			{
				std::cout << "# - WARNING: Unrecognized option \"" << strargv << "\" ignored!" << std::endl;
			}
		}
	}
	else
	{
		std::cout << "# - No commandline options specified." << std::endl;
	}
	// -----------------------------------------------------
	// END of commandline parsing
	// -----------------------------------------------------

	RunSection::RunSection rs;
	rs.SetOverruleAppend(appendMode);
	rs.SetNoCalculationsMode(noCalculations);
	// rs.SetThreads(threads);

	arma::wall_clock timer;
	timer.tic();

	// Load input file to setup the RunSection object
	{
		if (!silentMode)
		{
			std::cout << hline << std::endl;
			std::cout << "# Loading input file(s)..." << std::endl;
		}

		std::string strargv(argv[argc - 1]);
		MSDParser::MSDParser parser(strargv);

		// Attempt to load the input file
		if (!parser.Load())
		{
			std::cout << "ERROR: Failed to open file \"" << strargv << "\"!" << std::endl;
			return 1;
		}

		parser.FillRunSection(rs);
	}

	// Notify that we are ready to perform actual calculations
	if (!silentMode)
	{
		std::cout << hline << std::endl;
		std::cout << "# Input loaded and objects prepared in " << timer.toc() << " seconds." << std::endl;
		std::cout << hline << std::endl;
	}

	// Output all defined names and included files
	if (showDefines)
	{
		auto fl = MSDParser::FileReader::GetFileList();
		auto defnames = MSDParser::FileReader::GetDefinitions();

		std::cout << "# Included files:" << std::endl;
		for (auto i = fl.cbegin(); i != fl.cend(); i++)
			std::cout << " - " << (*i) << std::endl;

		std::cout << "\n# Defined names:" << std::endl;
		for (auto i = defnames.cbegin(); i != defnames.cend(); i++)
		{
			std::cout << " - " << i->first << ": ";
			if (i->second.empty())
				std::cout << "<empty>";
			else
				std::cout << i->second;
			std::cout << std::endl;
		}
		std::cout << hline << std::endl;
	}

	// Output all the objects that were read from the input file
	if (showObjects)
	{
		rs.PrintSystems(showFullObjects);
		std::cout << hline << std::endl;
	}

	// If a list of ActionTargets has been requested, print them
	if (printActionTargets)
	{
		auto actionScalars = rs.GetActionScalars();
		auto actionVectors = rs.GetActionVectors();
		std::cout << "# ActionScalars:" << std::endl;
		for (auto i = actionScalars.cbegin(); i != actionScalars.cend(); i++)
			std::cout << " - " << i->first << ((i->second.IsReadonly()) ? " (readonly)" : "") << std::endl;
		std::cout << "\n# ActionVectors:" << std::endl;
		for (auto i = actionVectors.cbegin(); i != actionVectors.cend(); i++)
			std::cout << " - " << i->first << ((i->second.IsReadonly()) ? " (readonly)" : "") << std::endl;
		std::cout << hline << std::endl;
	}

	// Definitions
	auto steps = rs.GetSettings()->Steps();
	double runtime = 0.0;
	double totalruntime = 0.0;

	// If we should do the calculations, do them
	if (!noCalculations)
	{
		// Perform steps if we should start later than at step 1
		for (unsigned int i = 1; i < firstStep && i <= steps; i++)
		{
			rs.Step(i + 1);
		}

		for (unsigned int i = firstStep; i <= steps; i++)
		{
			// Check that we have not exceeded the step limit (if any)
			if (stepLimit > 0 && i - firstStep >= stepLimit)
			{
				steps = stepLimit; // This is the number of steps that was run, and which is used to calculate the average time per step
				break;
			}

			// Information about the step we are about to run
			if (!silentMode && i % reportSteps == 0)
			{
				std::cout << "# Now running step " << i << "/" << steps << "." << std::endl;
				std::cout << hline << std::endl;
			}

			// Start the timer
			timer.tic();

			// Run all the tasks in the RunSection - skip some if we have a checkpoint
			if (hasCheckpoint)
			{
				rs.Run(checkpoint, i);
				hasCheckpoint = false; // We only use the checkpoint for one step
			}
			else
			{
				rs.Run(i);
			}

			// Advance to the next calculation step
			rs.Step(i + 1);

			// Get the time for the step
			runtime = timer.toc();
			totalruntime += runtime;

			// Show time for the step after it finished
			if (!silentMode && i % reportSteps == 0)
			{
				std::cout << hline << std::endl;
				std::cout << "# Finished with step " << i << "/" << steps << " in " << runtime << " seconds." << std::endl;
			}
		}

		std::cout << hline << std::endl;
		std::cout << "# Calculations done in " << totalruntime << " seconds, with an average runtime per step of " << (totalruntime / (double)steps) << " seconds." << std::endl;
	}
	else
	{
		std::cout << hline << std::endl;
		std::cout << "# Shutting down without doing calculations." << std::endl;
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
