/////////////////////////////////////////////////////////////////////////
// TaskStaticSSSpectra implementation (RunSection module)  developed by Irina Anisimova and Luca Gerhards.
//
// Molecular Spin Dynamics Software - developed by Luca Gerhards.
// (c) 2022 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include "TaskStaticSSSpectra.h"
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "Spin.h"
#include "Interaction.h"
#include "ObjectParser.h"
#include "Operator.h"
#include "Pulse.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskStaticSSSpectra Constructors and Destructor
	// -----------------------------------------------------
	TaskStaticSSSpectra::TaskStaticSSSpectra(const MSDParser::ObjectParser &_parser, const RunSection &_runsection) : BasicTask(_parser, _runsection), timestep(1.0), totaltime(1.0e+4), reactionOperators(SpinAPI::ReactionOperatorType::Haberkorn)
	{
	}

	TaskStaticSSSpectra::~TaskStaticSSSpectra()
	{
	}
	// -----------------------------------------------------
	// TaskStaticSSSpectra protected methods
	// -----------------------------------------------------
	bool TaskStaticSSSpectra::RunLocal()
	{
		this->Log() << "Running method StaticSS-Spectra." << std::endl;

		// If this is the first step, write first part of header to the data file
		if (this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}

		// Decline density matrix and density vector variables
		arma::cx_mat rho0;
		arma::cx_vec rho0vec;

		// Obtain spin systems
		auto systems = this->SpinSystems();

		// Loop through all SpinSystems
		for (auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Make sure we have an initial state
			auto initial_states = (*i)->InitialState();
			if (initial_states.size() < 1)
			{
				this->Log() << "Skipping SpinSystem \"" << (*i)->Name() << "\" as no initial state was specified." << std::endl;
				continue;
			}

			this->Log() << "\nStarting with SpinSystem \"" << (*i)->Name() << "\"." << std::endl;

			// Obtain a SpinSpace to describe the system
			SpinAPI::SpinSpace space(*(*i));
			space.UseSuperoperatorSpace(true);
			space.SetReactionOperatorType(this->reactionOperators);

			std::vector<double> weights;
			weights = (*i)->Weights();

			// Normalize the weights
			double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
			if (sum_weights > 0)
			{
				for (double &weight : weights)
				{
					weight /= sum_weights;
				}
			}

			// Get the initial state
			if (weights.size() > 1)
			{
				this->Log() << "Using weighted density matrix for initial state. Be sure that the sum of weights equals to 1." << std::endl;
				// Get the initial state
				int counter = 0;
				for (auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
				{
					arma::cx_mat tmp_rho0;
					if (!space.GetState(*j, tmp_rho0))
					{
						this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
						continue;
					}

					if (j == initial_states.cbegin())
					{
						this->Log() << "State: \"" << (*j)->Name() << "\", Weight:\"" << weights[0] << "\"." << std::endl;
						rho0 = weights[0] * tmp_rho0;
						counter += 1;
					}
					else
					{
						this->Log() << "State: \"" << (*j)->Name() << "\", Weight:\"" << weights[counter] << "\"." << std::endl;
						rho0 += weights[counter] * tmp_rho0;
						counter += 1;
					}
				}
			}
			else
			// Get the initial state without weights
			{
				for (auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
				{
					arma::cx_mat tmp_rho0;

					// Get the initial state in thermal equilibrium
					if ((*j) == nullptr) // "Thermal initial state"
					{
						this->Log() << "Initial state = thermal " << std::endl;

						// Get the thermalhamiltonianlist
						std::vector<std::string> thermalhamiltonian_list = (*i)->ThermalHamiltonianList();

						this->Log() << "ThermalHamiltonianList = [";
						for (size_t j = 0; j < thermalhamiltonian_list.size(); j++)
						{
							this->Log() << thermalhamiltonian_list[j];
							if (j < thermalhamiltonian_list.size() - 1)
								this->Log() << ", "; // Add a comma between elements
						}
						this->Log() << "]" << std::endl;

						// Get temperature
						double temperature = (*i)->Temperature();
						this->Log() << "Temperature = " << temperature << "K" << std::endl;

						// Get the initial state with thermal equilibrium
						if (!space.GetThermalState(space, temperature, thermalhamiltonian_list, tmp_rho0))
						{
							this->Log() << "Failed to obtain projection matrix onto thermal state, initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
							continue;
						}
					}
					else // Get the initial state without thermal equilibrium
					{
						if (!space.GetState(*j, tmp_rho0))
						{
							this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
							continue;
						}
					}

					// Obtain the initial density matrix
					if (j == initial_states.cbegin())
					{
						rho0 = tmp_rho0;
					}
					else
					{
						rho0 += tmp_rho0;
					}
				}
			}

			rho0 /= arma::trace(rho0); // The density operator should have a trace of 1

			// Convert initial state to superoperator space
			if (!space.OperatorToSuperspace(rho0, rho0vec))
			{
				this->Log() << "Failed to convert initial state density operator to superspace." << std::endl;
				continue;
			}

			// Get the Hamiltonian
			arma::sp_cx_mat H;
			if (!space.Hamiltonian(H))
			{
				this->Log() << "Failed to obtain Hamiltonian in superspace." << std::endl;
				continue;
			}

			// Get a matrix to collect all the terms (the total Liouvillian)
			arma::sp_cx_mat A = arma::cx_double(0.0, -1.0) * H;

			// Get the reaction operators, and add them to "A"
			arma::sp_cx_mat K;
			if (!space.TotalReactionOperator(K))
			{
				this->Log() << "Warning: Failed to obtain matrix representation of the reaction operators!" << std::endl;
			}
			A -= K;

			// Get the relaxation terms, assuming that they can just be added to the Liouvillian superoperator
			arma::sp_cx_mat R;
			for (auto j = (*i)->operators_cbegin(); j != (*i)->operators_cend(); j++)
			{
				if (space.RelaxationOperator((*j), R))
				{
					A += R;
					this->Log() << "Added relaxation operator \"" << (*j)->Name() << "\" to the Liouvillian.\n";
				}
			}

			// Read in input parameters 

			// Read the method from the input file
			std::string Method;
			if (!this->Properties()->Get("method", Method))
			{
				this->Log() << "Failed to obtain an input for a Method. Please specify method = timeinf or method = timeevo." << std::endl;
			}

			// Read if the result should be integrated or not if method.
			bool integration = false;
			if (!this->Properties()->Get("integration", integration))
			{
				this->Log() << "Failed to obtain an input for an integtation. Plese use integration = true/false. Using integration = false by default. " << std::endl;
			}
			this->Log() << "Integration of the yield in time on a grid  = " << integration << std::endl;

			// Read integrationwindow from the input file
			std::string Integrationwindow;
			if (!this->Properties()->Get("integrationtimeframe", Integrationwindow))
			{
				this->Log() << "Failed to obtain an input for a integrationtimeframe. Please choose integrationtimeframe =  pulse / freeevo / full. Using freeevo propagation evolution window by default" << std::endl;
				Integrationwindow = "freeevo";
			}
			this->Log() << "Timewindow for the propagation integration: " << Integrationwindow << std::endl;

			// Read CIDSP from the input file
			bool CIDSP = false;
			if (!this->Properties()->Get("cidsp", CIDSP))
			{
				this->Log() << "Failed to obtain an input for a CIDSP. Plese use cidsp = true/false. Using cidsp = false by default. " << std::endl;
			}

			// Read printtimeframe from the input file
			std::string Timewindow;
			if (!this->Properties()->Get("printtimeframe", Timewindow))
			{
				this->Log() << "Failed to obtain an input for a printtimeframe. Please choose printtimeframe =  pulse / freeevo / full. Using full propagation evolution window by default" << std::endl;\
				Timewindow = "full";
			}
			this->Log() << "Timewindow for the propagation printing: " << Timewindow << std::endl;

			double Printedtime = 0;

			// Initialize a first step
			arma::cx_vec rhovec = rho0vec;

			// Read a pulse sequence from the input
			std::vector<std::tuple<std::string, double>> Pulsesequence;
			if (this->Properties()->GetPulseSequence("pulsesequence", Pulsesequence))
			{
				this->Log() << "Pulsesequence" << std::endl;

				// Loop through all pulse sequences
				for (const auto &seq : Pulsesequence)
				{
					// Write which pulse in pulsesequence is calculating now
					this->Log() << std::get<0>(seq) << ", " << std::get<1>(seq) << std::endl;

					// Save the parameters from the input as variables
					std::string pulse_name = std::get<0>(seq);
					double timerelaxation = std::get<1>(seq);

					for (auto pulse = (*i)->pulses_cbegin(); pulse < (*i)->pulses_cend(); pulse++)
					{
						if ((*pulse)->Name().compare(pulse_name) == 0)
						{
							// Apply a pulse to our density vector
							if ((*pulse)->Type() == SpinAPI::PulseType::InstantPulse)
							{
								// Create a Pulse operator in SS
								arma::sp_cx_mat pulse_operator;
								if (!space.PulseOperator((*pulse), pulse_operator))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}
								rhovec = pulse_operator * rhovec;
							}
							else if ((*pulse)->Type() == SpinAPI::PulseType::LongPulseStaticField)
							{
								// Create a Pulse operator in SS
								arma::sp_cx_mat pulse_operator;
								if (!space.PulseOperator((*pulse), pulse_operator))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}

								// Create a holder vector for an averaged density
								arma::cx_vec rhoavg;
								rhoavg.zeros(size(rho0vec));

								int firststep;
								if (Printedtime == 0)
									firststep = 0;
								else
									firststep = 1;

								// Create array containing a propagator and the current state of each system
								std::pair<arma::sp_cx_mat, arma::cx_vec> G;
								// Get the propagator and put it into the array together with the initial state
								arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from((A + (arma::cx_double(0.0, -1.0) * pulse_operator)) * (*pulse)->Timestep())));
								G = std::pair<arma::sp_cx_mat, arma::cx_vec>(A_sp, rhovec);

								unsigned int steps = static_cast<unsigned int>(std::abs((*pulse)->Pulsetime() / (*pulse)->Timestep()));
								for (unsigned int n = firststep; n <= steps; n++)
								{
									if (!n == 0) // for the very first step just printing is done, no propagation
									{
										// Take a step, "first" is propagator and "second" is current state
										rhovec = G.first * G.second;

										// Integrate the result if needed
										if ((integration) && (Integrationwindow.compare("freeevo") != 0))
										{
											rhoavg += (*pulse)->Timestep() * (G.second + rhovec) / 2;
										}

										// Get the new current state density vector
										G.second = rhovec;

										// Save the result if there were some changes
										if (!rhoavg.is_zero(0))
										{
											rhovec = rhoavg;
										}
									}

									if (Timewindow.compare("freeevo") != 0)
									{
										if (!this->ProjectAndPrintOutputLine(i, space, rhovec, Printedtime, (*pulse)->Timestep(), n, CIDSP, this->Data(), this->Log()))
											this->Log() << "Could not project the state vector and print the result into a file" << std::endl;
									}
								}
							}
							else if ((*pulse)->Type() == SpinAPI::PulseType::LongPulse)
							{
								// Create a Pulse operator in SS
								arma::sp_cx_mat pulse_operator;
								if (!space.PulseOperator((*pulse), pulse_operator))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}

								// Create a holder vector for an averaged density
								arma::cx_vec rhoavg;
								rhoavg.zeros(size(rho0vec));

								int firststep;
								if (Printedtime == 0)
									firststep = 0;
								else
									firststep = 1;

								// Create array containing a propagator and the current state of each system
								std::pair<arma::sp_cx_mat, arma::cx_vec> G;

								// Get the propagator and put it into the array together with the initial state
								arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from((A + (arma::cx_double(0.0, -1.0) * pulse_operator * std::cos((*pulse)->Frequency() * (*pulse)->Timestep()))) * (*pulse)->Timestep())));
								G = std::pair<arma::sp_cx_mat, arma::cx_vec>(A_sp, rhovec);

								unsigned int steps = static_cast<unsigned int>(std::abs((*pulse)->Pulsetime() / (*pulse)->Timestep()));
								for (unsigned int n = firststep; n <= steps; n++)
								{
									if (!n == 0) // for the very first step just printing is done, no propagation
									{
										// Take a step, "first" is propagator and "second" is current state
										rhovec = G.first * G.second;

										// Integrate the result if needed
										if ((integration) && (Integrationwindow.compare("freeevo") != 0))
										{
											rhoavg += (*pulse)->Timestep() * (G.second + rhovec) / 2;
										}

										// Get the new current state density vector
										G.second = rhovec;

										// Save the result if there were some changes
										if (!rhoavg.is_zero(0))
										{
											rhovec = rhoavg;
										}
									}

									if (Timewindow.compare("freeevo") != 0)
									{
										if (!this->ProjectAndPrintOutputLine(i, space, rhovec, Printedtime, (*pulse)->Timestep(), n, CIDSP, this->Data(), this->Log()))
											this->Log() << "Could not project the state vector and print the result into a file" << std::endl;
									}

								}
							}
							else if ((*pulse)->Type() == SpinAPI::PulseType::MWPulse)
							{
								// Create a Pulse operator in SS
								arma::sp_cx_mat pulse_operator;
								arma::sp_cx_mat lrot;
								arma::sp_cx_mat lroths;
								if (!space.PulseOperator_mw((*pulse), pulse_operator, lrot, lroths))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}

								// Create a holder vector for an averaged density
								arma::cx_vec rhoavg;
								rhoavg.zeros(size(rho0vec));

								int firststep;
								if (Printedtime == 0)
									firststep = 0;
								else
									firststep = 1;

								// Create array containing a propagator and the current state of each system
								std::pair<arma::sp_cx_mat, arma::cx_vec> G;

								// Get the propagator and put it into the array together with the initial state
								arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from((A * (*pulse)->Timestep()) + (lrot * pulse_operator * lrot.t()))));
								// arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from((A * (*pulse)->Timestep()) + pulse_operator)));
								G = std::pair<arma::sp_cx_mat, arma::cx_vec>(A_sp, rhovec);

								unsigned int steps = static_cast<unsigned int>(std::abs((*pulse)->Pulsetime() / (*pulse)->Timestep()));
								for (unsigned int n = firststep; n <= steps; n++)
								{
									if (!n == 0) // for the very first step just printing is done, no propagation
									{
										// Take a step, "first" is propagator and "second" is current state
										rhovec = G.first * G.second;

										// Integrate the result if needed
										if ((integration) && (Integrationwindow.compare("freeevo") != 0))
										{
											rhoavg += (*pulse)->Timestep() * (G.second + rhovec) / 2;
										}

										// Get the new current state density vector
										G.second = rhovec;

										// Save the result if there were some changes
										if (!rhoavg.is_zero(0))
										{
											rhovec = rhoavg;
										}
									}

									// RK45
									// // arma::sp_cx_mat A_sp = (A*(*pulse)->Timestep()) + pulse_operator;
									// arma::sp_cx_mat A_sp = lrot * ((A * (*pulse)->Timestep()) + pulse_operator) * lrot.t();
									// // rhovec = lrot * rhovec;

									// double InitialTimestep = (*pulse)->Timestep();
									// double Currenttime = 0;
									// double time = 0;
									// unsigned int n = 0;
									// while (Currenttime <= (*pulse)->Pulsetime())
									// {
									// 	// std::cout << n << Currenttime << std::endl;

									// 	if (!n == 0)
									// 	{
									// 		{
									// 			Currenttime += (*pulse)->Timestep();
									// 			time = RungeKutta45Armadillo(A_sp, rhovec, rhovec, (*pulse)->Timestep(), ComputeRhoDot, {1e-7, 1e-6}, InitialTimestep * 1e-3, InitialTimestep * 1e4);
									// 		}
									// 	}

									// // rhovec = lrot.t() * rhovec

									if (Timewindow.compare("freeevo") != 0)
									{
										if (!this->ProjectAndPrintOutputLine(i, space, rhovec, Printedtime, (*pulse)->Timestep(), n, CIDSP, this->Data(), this->Log()))
											this->Log() << "Could not project the state vector and print the result into a file" << std::endl;
									}

									// n++; //RK45
								}
							}
							else
							{
								this->Log() << "Current pulse type is not implemented. Please use type = InstantPulse / LongPulseStaticField / LongPulse / MWPulse." << std::endl;
							}

							// Update the printed time according to the printtimeframe key
							if (Timewindow.compare("freeevo") != 0)
							{
								Printedtime += (*pulse)->Pulsetime();
							}

							// Get the system relax during the time

							if (!timerelaxation == 0)
							{
								// Create a holder vector for an averaged density
								arma::cx_vec rhoavg;
								rhoavg.zeros(size(rho0vec));

								// Create array containing a propagator and the current state of each system
								std::pair<arma::sp_cx_mat, arma::cx_vec> G;
								arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from(A * (*pulse)->Timestep())));
								// Get the propagator and put it into the array together with the initial state
								G = std::pair<arma::sp_cx_mat, arma::cx_vec>(A_sp, rhovec);

								unsigned int steps = static_cast<unsigned int>(std::abs(timerelaxation / (*pulse)->Timestep()));
								for (unsigned int n = 1; n <= steps; n++) // n always starts with one here, because this part cannot be done without the pulse part on top
								{
									// Take a step, "first" is propagator and "second" is current state
									rhovec = G.first * G.second;

									// Integrate the result if needed
									if ((integration) && (Integrationwindow.compare("freeevo") != 0))
									{
										rhoavg += (*pulse)->Timestep() * (G.second + rhovec) / 2;
									}

									// Get the new current state density vector
									G.second = rhovec;

									// Save the result if there were some changes
									if (!rhoavg.is_zero(0))
									{
										rhovec = rhoavg;
									}

									if (Timewindow.compare("freeevo") != 0)
									{
										if (!this->ProjectAndPrintOutputLine(i, space, rhovec, Printedtime, (*pulse)->Timestep(), n, CIDSP, this->Data(), this->Log()))
											this->Log() << "Could not project the state vector and print the result into a file" << std::endl;
									}
								}
							}
							// Update the printed time according to the printtimeframe key
							if (Timewindow.compare("freeevo") != 0)
							{
								Printedtime += timerelaxation;
							}
						}
					}
				}
			}

			// Free time evolution propagation after the pulse sequience

			// Method Propagation to infinity
			if (Method.compare("timeinf") == 0)
			{
				// Perform the calculation
				this->Log() << "Ready to perform calculation." << std::endl;

				this->Log() << "Method = " << Method << std::endl;

				// Print the warning to the user, what integration keyword is used
				if (integration)
				{
					this->Log() << "Warning: steady state method (timeinf) is calculated as an inverse of the Liouvillian operator, instead of the integration on a grid."<<
					"The integration of the pulse sequence timewindow could be added if integration = true and integrationtimeframe = pulse / full." << std::endl;
				}

				arma::cx_vec result = -solve(arma::conv_to<arma::cx_mat>::from(A), rhovec);

				rhovec = result;

				if (Timewindow.compare("pulse") != 0)
				{
					if (!this->ProjectAndPrintOutputLineInf(i, space, rhovec, Printedtime, this->timestep, CIDSP, this->Data(), this->Log()))
						this->Log() << "Could not project the state vector and print the result into a file" << std::endl;
				}

				this->Log() << "Done with calculation." << std::endl;
			}
			// Method TIME EVOLUTION
			else if (Method.compare("timeevo") == 0)
			{

				if (!this->totaltime == 0)
				{
					// Perform the calculation
					this->Log() << "Ready to perform calculation." << std::endl;

					this->Log() << "Method = " << Method << std::endl;

					// Create a holder vector for an averaged density
					arma::cx_vec rhoavg;
					rhoavg.zeros(size(rho0vec));

					// Avoid printing double timesteps
					int firststep;
					if (Printedtime == 0)
						firststep = 0;
					else
						firststep = 1;

					// Create array containing a propagator and the current state of each system
					std::pair<arma::sp_cx_mat, arma::cx_vec> G;
					arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from(A * this->timestep)));
					// Get the propagator and put it into the array together with the initial state
					G = std::pair<arma::sp_cx_mat, arma::cx_vec>(A_sp, rhovec);

					unsigned int steps = static_cast<unsigned int>(std::abs(this->totaltime / this->timestep));
					for (unsigned int n = firststep; n <= steps; n++)
					{
						if (!n == 0)
						{
							// Take a step, "first" is propagator and "second" is current state
							rhovec = G.first * G.second;

							// Integrate the density vector over the current time interval
							if ((integration) && (Integrationwindow.compare("pulse") != 0))
							{
								rhoavg += this->timestep * (G.second + rhovec) / 2;
							}

							// Get the new current state density vector
							G.second = rhovec;

							// Save the result if there were some changes
							if (!rhoavg.is_zero(0))
							{
								rhovec = rhoavg;
							}
						}

						// double InitialTimestep = this->timestep;
						// double Currenttime = 0;
						// double time;
						// unsigned int n = 0;
						// while (Currenttime <= this->totaltime)
						// {
						// 	if (!n == 0)
						// 	{
						// 		{
						// 			Currenttime += this->timestep;
						// 			time = RungeKutta45Armadillo(A, rhovec, rhovec, this->timestep, ComputeRhoDot, {1e-7, 1e-6}, InitialTimestep * 1e-3, InitialTimestep * 1e4);
						// 		}
						// 	}

						if (Timewindow.compare("pulse") != 0)
						{
							if (!this->ProjectAndPrintOutputLine(i, space, rhovec, Printedtime, this->timestep, n, CIDSP, this->Data(), this->Log()))
								this->Log() << "Could not project the state vector and print the result into a file" << std::endl;
						}

						// n++; //delete if not RK
					}
				}

				this->Log() << "Done with calculation." << std::endl;
			}
			else
			{
				this->Log() << "Undefined spectroscopy method. Please choose between timeinf or timeevo methods." << std::endl;
			}

			this->Log() << "\nDone with SpinSystem \"" << (*i)->Name() << "\"" << std::endl;
		}

		// Terminate the line in the data file after iteration through all spin systems
		this->Data() << std::endl;

		return true;
	}

	// Writes the header of the data file (but can also be passed to other streams)
	void TaskStaticSSSpectra::WriteHeader(std::ostream &_stream)
	{
		_stream << "Step ";
		_stream << "Time ";
		this->WriteStandardOutputHeader(_stream);

		std::vector<std::string> spinList;
		bool CIDSP = false;
		int m;

		// Get header for each spin system
		auto systems = this->SpinSystems();
		for (auto i = systems.cbegin(); i != systems.cend(); i++)
		{

			if (this->Properties()->GetList("spinlist", spinList, ','))
			{
				for (auto l = (*i)->spins_cbegin(); l != (*i)->spins_cend(); l++)
				{
					std::string spintype;

					(*l)->Properties()->Get("type", spintype);

					for (m = 0; m < (int)spinList.size(); m++)
					{

						if ((*l)->Name() == spinList[m])
						{
							// Yields are written per transition
							// bool CIDSP = false;
							if (this->Properties()->Get("cidsp", CIDSP) && CIDSP == true)
							{
								// Write each transition name
								auto transitions = (*i)->Transitions();
								for (auto j = transitions.cbegin(); j != transitions.cend(); j++)
								{
									_stream << (*i)->Name() << "." << (*l)->Name() << "." << (*j)->Name() << ".yield"
											<< ".Ix ";
									_stream << (*i)->Name() << "." << (*l)->Name() << "." << (*j)->Name() << ".yield"
											<< ".Iy ";
									_stream << (*i)->Name() << "." << (*l)->Name() << "." << (*j)->Name() << ".yield"
											<< ".Iz ";
								}
							}
							else
							{
								// Write each state name
								auto states = (*i)->States();
								_stream << (*i)->Name() << "." << (*l)->Name() << ".Ix ";
								_stream << (*i)->Name() << "." << (*l)->Name() << ".Iy ";
								_stream << (*i)->Name() << "." << (*l)->Name() << ".Iz ";
							}
						}
					}
				}
			}
		}
		_stream << std::endl;
	}

	// Validation
	bool TaskStaticSSSpectra::Validate()
	{

		double inputTimestep = 0.0;
		double inputTotaltime = 0.0;

		// Get timestep
		if (this->Properties()->Get("timestep", inputTimestep))
		{
			if (std::isfinite(inputTimestep) && inputTimestep > 0.0)
			{
				this->timestep = inputTimestep;
			}
			else
			{
				this->Log() << "# WARNING: undefined timestep, using by default 0.1 ns!" << std::endl;
				this->timestep = 0.1;
			}
		}

		// Get totaltime
		if (this->Properties()->Get("totaltime", inputTotaltime))
		{
			if (std::isfinite(inputTotaltime) && inputTotaltime >= 0.0)
			{
				this->totaltime = inputTotaltime;
			}
			else
			{
				this->Log() << "# ERROR: invalid total time!" << std::endl;
				return false;
			}
		}

		// Get the reacton operator type
		std::string str;
		if (this->Properties()->Get("reactionoperators", str))
		{
			if (str.compare("haberkorn") == 0)
			{
				this->reactionOperators = SpinAPI::ReactionOperatorType::Haberkorn;
				this->Log() << "Setting reaction operator type to Haberkorn." << std::endl;
			}
			else if (str.compare("lindblad") == 0)
			{
				this->reactionOperators = SpinAPI::ReactionOperatorType::Lindblad;
				this->Log() << "Setting reaction operator type to Lindblad." << std::endl;
			}
			else
			{
				this->Log() << "Warning: Unknown reaction operator type specified. Using default reaction operators." << std::endl;
			}
		}

		return true;
	}

	bool TaskStaticSSSpectra::GetEigenvectors_H0(SpinAPI::SpinSpace &_space, arma::vec &_eigen_val, arma::sp_cx_mat &_eigen_vec_sp) const
	{
		_space.UseSuperoperatorSpace(false);

		arma::cx_mat H;

		if (!_space.Hamiltonian(H))
		{
			// this->Log() << "Failed to obtain Static Hamiltonian in Hilbert Space." << std::endl;
		}

		arma::cx_mat _eigen_vec;

		// ----------------------------------------------------------------
		// DIAGONALIZATION OF H0// We need all of these operators
		// ----------------------------------------------------------------
		// this->Log() << "Starting diagonalization..." << std::endl;
		arma::eig_sym(_eigen_val, _eigen_vec, (H));
		// this->Log() << "Diagonalization done! Eigenvalues: " << _eigen_val.n_elem << ", eigenvectors: " << _eigen_vec.n_cols << std::endl;

		_eigen_vec_sp = arma::conv_to<arma::sp_cx_mat>::from(_eigen_vec);

		_space.UseSuperoperatorSpace(true);

		return true;
	}

	arma::cx_vec TaskStaticSSSpectra::ComputeRhoDot(double t, arma::sp_cx_mat &L, arma::cx_vec &K, arma::cx_vec RhoNaught)
	{
		arma::cx_vec ReturnVec(L.n_rows);
		RhoNaught = RhoNaught + K;
		ReturnVec = L * RhoNaught;
		return ReturnVec;
	}

	bool TaskStaticSSSpectra::ProjectAndPrintOutputLine(auto &_i, SpinAPI::SpinSpace &_space, arma::cx_vec &_rhovec, double &_printedtime, double _timestep, unsigned int &_n, bool &_cidsp, std::ostream &_datastream, std::ostream &_logstream)
	{
		arma::cx_mat rho0;

		// Convert the resulting density operator back to its Hilbert space representation
		if ((!_space.OperatorFromSuperspace(_rhovec, rho0)) && (_n == 0))
		{
			_logstream << "Failed to convert resulting superspace-vector back to native Hilbert space." << std::endl;
			return false;
		}

		// Get nuclei of interest for CIDNP spectrum
		arma::cx_mat Iprojx;
		arma::cx_mat Iprojy;
		arma::cx_mat Iprojz;

		std::vector<std::string> spinList;

		if (_n == 0)
			_logstream << "CIDSP = " << _cidsp << std::endl;

		// Save the current step
		_datastream << this->RunSettings()->CurrentStep() << " ";
		// Save the current time
		_datastream << std::setprecision(12) << _printedtime + (_n * _timestep) << " ";
		this->WriteStandardOutput(_datastream);

		if (this->Properties()->GetList("spinlist", spinList, ','))
		{

			for (auto l = (*_i)->spins_cbegin(); l != (*_i)->spins_cend(); l++)
			{
				for (int m = 0; m < (int)spinList.size(); m++)
				{
					if ((*l)->Name() == spinList[m])
					{
						if (!_space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*l)->Sx()), (*l), Iprojx))
						{
							return false;
						}

						if (!_space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*l)->Sy()), (*l), Iprojy))
						{
							return false;
						}

						if (!_space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*l)->Sz()), (*l), Iprojz))
						{
							return false;
						}


						arma::cx_mat P;

						// There are two result modes - either write results per transition  if CIDSP is true or for each defined state if CIDSP is false

						if (_cidsp == true)
						{
							// Loop through all defind transitions
							auto transitions = (*_i)->Transitions();
							for (auto j = transitions.cbegin(); j != transitions.cend(); j++)
							{
								// Make sure that there is a state object
								if ((*j)->SourceState() == nullptr)
									continue;

								if ((!_space.GetState((*j)->SourceState(), P)) && (_n == 0))
								{
									_logstream << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*_i)->Name() << "\"." << std::endl;
									continue;
								}


								// Return the yield for this transition
								_datastream << std::real(arma::trace(Iprojx * (*j)->Rate() * P * rho0)) << " ";
								_datastream << std::real(arma::trace(Iprojy * (*j)->Rate() * P * rho0)) << " ";
								_datastream << std::real(arma::trace(Iprojz * (*j)->Rate() * P * rho0)) << " ";
							}
						}
						else if (_cidsp == false)
						{
							// Return the yield for this state - note that no reaction rates are included here.
							_datastream << std::real(arma::trace(Iprojx * rho0)) << " ";
							_datastream << std::real(arma::trace(Iprojy * rho0)) << " ";
							_datastream << std::real(arma::trace(Iprojz * rho0)) << " ";
						}
					}
				}
			}
		}
		else
		{
			if (_n == 0)
				_logstream << "No nucleus was specified for projection" << std::endl;
			return false;
		}

		_datastream << std::endl;

		return true;
	}

	bool TaskStaticSSSpectra::ProjectAndPrintOutputLineInf(auto &_i, SpinAPI::SpinSpace &_space, arma::cx_vec &_rhovec, double &_printedtime, double _timestep, bool &_cidsp, std::ostream &_datastream, std::ostream &_logstream)
	{
		arma::cx_mat rho0;

		// Convert the resulting density operator back to its Hilbert space representation
		if ((!_space.OperatorFromSuperspace(_rhovec, rho0)))
		{
			_logstream << "Failed to convert resulting superspace-vector back to native Hilbert space." << std::endl;
			return false;
		}

		// Get nuclei of interest for CIDNP spectrum
		arma::cx_mat Iprojx;
		arma::cx_mat Iprojy;
		arma::cx_mat Iprojz;

		std::vector<std::string> spinList;

		_logstream << "CIDSP = " << _cidsp << std::endl;

		// Save the current step
		_datastream << this->RunSettings()->CurrentStep() << " ";
		// Save the current time
		_datastream << "inf" << " ";
		this->WriteStandardOutput(_datastream);

		if (this->Properties()->GetList("spinlist", spinList, ','))
		{

			for (auto l = (*_i)->spins_cbegin(); l != (*_i)->spins_cend(); l++)
			{
				for (int m = 0; m < (int)spinList.size(); m++)
				{
					if ((*l)->Name() == spinList[m])
					{
						if (!_space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*l)->Sx()), (*l), Iprojx))
						{
							return false;
						}

						if (!_space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*l)->Sy()), (*l), Iprojy))
						{
							return false;
						}

						if (!_space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*l)->Sz()), (*l), Iprojz))
						{
							return false;
						}

						arma::cx_mat P;

						// There are two result modes - either write results per transition  if CIDSP is true or for each defined state if CIDSP is false

						if (_cidsp == true)
						{
							// Loop through all defind transitions
							auto transitions = (*_i)->Transitions();
							for (auto j = transitions.cbegin(); j != transitions.cend(); j++)
							{
								// Make sure that there is a state object
								if ((*j)->SourceState() == nullptr)
									continue;

								if ((!_space.GetState((*j)->SourceState(), P)))
								{
									_logstream << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*_i)->Name() << "\"." << std::endl;
									continue;
								}

								// Return the yield for this transition
								_datastream << std::real(arma::trace(Iprojx * (*j)->Rate() * P * rho0)) << " ";
								_datastream << std::real(arma::trace(Iprojy * (*j)->Rate() * P * rho0)) << " ";
								_datastream << std::real(arma::trace(Iprojz * (*j)->Rate() * P * rho0)) << " ";
							}
						}
						else if (_cidsp == false)
						{
							// Return the yield for this state - note that no reaction rates are included here.
							_datastream << std::real(arma::trace(Iprojx * rho0)) << " ";
							_datastream << std::real(arma::trace(Iprojy * rho0)) << " ";
							_datastream << std::real(arma::trace(Iprojz * rho0)) << " ";
						}
					}
				}
			}
		}
		else
		{
			_logstream << "No nucleus was specified for projection" << std::endl;
			return false;
		}

		return true;
	}

	// -----------------------------------------------------
}
