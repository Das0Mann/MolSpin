/////////////////////////////////////////////////////////////////////////
// TaskStaticHSDirectSpectra implementation (RunSection module) by Luca Gerhards
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TaskStaticHSDirectSpectra.h"
#include "Transition.h"
#include "Operator.h"
#include "Settings.h"
#include "State.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "ObjectParser.h"
#include "Spin.h"
#include "Interaction.h"
#include "Pulse.h"
#include <iomanip> // std::setprecision

namespace RunSection
{
	// -----------------------------------------------------
	// TaskStaticHSDirectSpectra Constructors and Destructor
	// -----------------------------------------------------
	TaskStaticHSDirectSpectra::TaskStaticHSDirectSpectra(const MSDParser::ObjectParser &_parser, const RunSection &_runsection) : BasicTask(_parser, _runsection), timestep(0.1), totaltime(1000), reactionOperators(SpinAPI::ReactionOperatorType::Haberkorn)
	{
	}

	TaskStaticHSDirectSpectra::~TaskStaticHSDirectSpectra()
	{
	}
	// -----------------------------------------------------
	// TaskStaticHSDirectSpectra protected methods
	// -----------------------------------------------------
	bool TaskStaticHSDirectSpectra::RunLocal()
	{
		this->Log() << "Running method StaticHS_Direct_Yields." << std::endl;

		// If this is the first step, write first part of header to the data file
		if (this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}

		// Loop through all SpinSystems
		auto systems = this->SpinSystems();
		for (auto i = systems.cbegin(); i != systems.cend(); i++) // iteration through all spin systems, in this case (or usually), this is one
		{
			// Count the number of nuclear spins
			int nucspins = 0;
			std::vector<int> SpinNumbers;
			for (auto l = (*i)->spins_cbegin(); l != (*i)->spins_cend(); l++)
			{
				std::string spintype;
				(*l)->Properties()->Get("type", spintype);
				if (spintype != "electron")
				{
					nucspins += 1;
				}
				if (spintype == "electron")
				{
					// Throws an error if the spins are not spin 1/2
					if ((*l)->Multiplicity() != 2)
					{
						std::cout << (*l)->Multiplicity() << std::endl;
						this->Log() << "SkippingSpin System \"" << (*i)->Name() << "\" as electron spins have the wrong multiplicity." << std::endl;
						std::cout << "# ERROR: electron spins have to be spin 1/2! Skipping the SpinSystem." << std::endl;
						return 1;
					}
				}
			}

			this->Log() << "\nStarting with SpinSystem \"" << (*i)->Name() << "\"." << std::endl;

			// Obtain a SpinSpace to describe the system
			SpinAPI::SpinSpace space(*(*i));
			space.UseSuperoperatorSpace(false);
			space.SetReactionOperatorType(this->reactionOperators);

			std::string InitialState;
			arma::cx_mat InitialStateVector;
			if(this->Properties()->Get("initialstate", InitialState))
			{
				// Set up states for time-propagation
				arma::cx_mat TaskInitialStateVector(4, 1);
				std::string InitialStateLower;

				// Convert the string to lowercase for case-insensitive comparison
				InitialStateLower.resize(InitialState.size());
				std::transform(InitialState.begin(), InitialState.end(), InitialStateLower.begin(), ::tolower);

				if (InitialStateLower == "singlet")
				{
					arma::cx_mat SingletState(4, 1);
					SingletState(0) = 0.0;
					SingletState(1) = 1.0 / sqrt(2);
					SingletState(2) = -1.0 / sqrt(2);
					SingletState(3) = 0.0;
					TaskInitialStateVector = SingletState;
					this->Log() << "Singlet initial state." << std::endl;
				}
				else if (InitialStateLower == "tripletminus")
				{
					arma::cx_mat TripletMinusState(4, 1);
					TripletMinusState(0) = 0.0;
					TripletMinusState(1) = 0.0;
					TripletMinusState(2) = 0.0;
					TripletMinusState(3) = 1.0;
					TaskInitialStateVector = TripletMinusState;
					this->Log() << "Triplet minus initial state." << std::endl;
				}
				else if (InitialStateLower == "tripletzero")
				{
					arma::cx_mat TripletZeroState(4, 1);
					TripletZeroState(0) = 0.0;
					TripletZeroState(1) = 1.0 / sqrt(2);
					TripletZeroState(2) = 1.0 / sqrt(2);
					TripletZeroState(3) = 0.0;
					TaskInitialStateVector = TripletZeroState;
					this->Log() << "Triplet zero initial state." << std::endl;
				}
				else if (InitialStateLower == "tripletplus")
				{
					arma::cx_mat TripletPlusState(4, 1);
					TripletPlusState(0) = 1.0;
					TripletPlusState(1) = 0.0;
					TripletPlusState(2) = 0.0;
					TripletPlusState(3) = 0.0;
					TaskInitialStateVector = TripletPlusState;
					this->Log() << "Triplet plus initial state." << std::endl;
				}
				else
				{
					std::cout << "# ERROR: Invalid initial state value! It is set to a Singlet state." << std::endl;
					this->Log() << "Initial state is undefined. Setting it to a Singlet state" << std::endl;
					arma::cx_mat SingletState(4, 1);
					SingletState(0) = 0.0;
					SingletState(1) = 1.0 / sqrt(2);
					SingletState(2) = -1.0 / sqrt(2);
					SingletState(3) = 0.0;
					TaskInitialStateVector = SingletState;
				}
				InitialStateVector = TaskInitialStateVector;
			}
			else
			{
				// Make sure we have an initial state
				auto initial_states = (*i)->InitialState();
				if (initial_states.size() < 1)
				{
					this->Log() << "Skipping SpinSystem \"" << (*i)->Name() << "\" as no initial state was specified." << std::endl;
					continue;
				}

				arma::cx_vec tmp_InitialStateVector;

				for (auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
				{
					if (!space.GetStateSubSpace(*j, tmp_InitialStateVector))
					{
						this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
						continue;
					}
				}

				InitialStateVector = arma::reshape(tmp_InitialStateVector, tmp_InitialStateVector.n_elem, 1);
			}

			int Z = space.SpaceDimensions() / InitialStateVector.n_rows; // Size of the nuclear spin subspace
			std::cout << "# Hilbert Space Size " << InitialStateVector.n_rows * Z << " x " << InitialStateVector.n_rows * Z << std::endl;
			this->Log() << "Hilbert Space Size " << InitialStateVector.n_rows * Z << " x " << InitialStateVector.n_rows * Z << std::endl;
			this->Log() << "Size of Nuclear Spin Subspace " << Z << std::endl;

			arma::cx_mat B;
			B.zeros(Z * InitialStateVector.n_rows, Z);

			for (int it = 0; it < Z; it++)
			{
				arma::colvec temp(Z);
				temp(it) = 1;
				B.col(it) = arma::kron(InitialStateVector, temp);
			}

			// Get Information about the polarization of choice
			bool CIDSP = false;
			if (!this->Properties()->Get("cidsp", CIDSP))
			{
				this->Log() << "Failed to obtain an input for a CIDSP" << std::endl;
			}

			// Get projectors of interest of the spectrum
			arma::sp_cx_mat Iprojx;
			arma::sp_cx_mat Iprojy;
			arma::sp_cx_mat Iprojz;

			std::vector<std::string> spinList;
			int m;

			// Check transitions, rates and projection operators
			auto transitions = (*i)->Transitions();
			arma::sp_cx_mat P;
			int num_transitions = 0;

			int projection_counter = 0;
			std::map<int, arma::sp_cx_mat> Operators;
			arma::vec rates(1, 1);

			// Getting the projection operators
			if (this->Properties()->GetList("spinlist", spinList, ','))
			{
				for (auto l = (*i)->spins_cbegin(); l != (*i)->spins_cend(); l++)
				{
					for (m = 0; m < (int)spinList.size(); m++)
					{
						if ((*l)->Name() == spinList[m])
						{
							if (!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from((*l)->Sx()), (*l), Iprojx))
							{
								return false;
							}
								if (!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from((*l)->Sy()), (*l), Iprojy))
							{
								return false;
							}
								if (!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from((*l)->Sz()), (*l), Iprojz))
							{
								return false;
							}

							if (CIDSP == true)
							{
								// Gather rates and operators
								for (auto j = transitions.cbegin(); j != transitions.cend(); j++)
								{
									if ((*j)->SourceState() == nullptr)
										continue;
									if (!space.GetState((*j)->SourceState(), P))
									{
										std::cout << "# ERROR: Could not obtain projection matrix!" << std::endl;
										this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
										return 1;
									}
									if (num_transitions != 0)
									{
										rates.insert_rows(num_transitions, 1);
									}

									Operators[projection_counter] = (*j)->Rate() * Iprojx * P;
									projection_counter += 1;
									Operators[projection_counter] = (*j)->Rate() * Iprojy * P;
									projection_counter += 1;
									Operators[projection_counter] = (*j)->Rate() * Iprojz * P;
									projection_counter += 1;

									rates(num_transitions) = (*j)->Rate();
									num_transitions++;
								}
							}
							else
							{
								// Gather rates and operators
								Operators[projection_counter] = Iprojx;
								projection_counter += 1;
								Operators[projection_counter] = Iprojy;
								projection_counter += 1;
								Operators[projection_counter] = Iprojz;
								projection_counter += 1;
							
								for (auto j = transitions.cbegin(); j != transitions.cend(); j++)
								{
									if ((*j)->SourceState() == nullptr)
										continue;

									if (num_transitions != 0)
									{
										rates.insert_rows(num_transitions, 1);
									}

									rates(num_transitions) = (*j)->Rate();
									num_transitions++;
								}
							}
						}	
					}
				}
			}

			// Get the Hamiltonian
			arma::sp_cx_mat H;
			if (!space.Hamiltonian(H))
			{
				this->Log() << "Failed to obtain the Hamiltonian in Hilbert Space." << std::endl;
				std::cout << "# ERROR: Failed to obtain the Hamiltonian!" << std::endl;
				return 1;
			}

			arma::sp_cx_mat K;
			K.zeros(InitialStateVector.n_rows * Z, InitialStateVector.n_rows * Z);

			for (auto j = transitions.cbegin(); j != transitions.cend(); j++)
			{
				if ((*j)->SourceState() == nullptr)
					continue;
				space.GetState((*j)->SourceState(), P);
				K += (*j)->Rate() / 2 * P;
			}

			// Get pulses and pulse the system
			arma::sp_cx_mat A = arma::cx_double(0.0, -1.0) * H - K;

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
								// Create a Pulse operator in HS; only one side of exponentials as we only propagate wavevectors
								arma::sp_cx_mat pulse_operator;
								if (!space.PulseOperatorOnStatevector((*pulse), pulse_operator))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}

								// Take a step, "first" is propagator and "second" is current state
								for (int i = 0; i < int(B.n_cols); i++) 
								{
									B.col(i) = pulse_operator * B.col(i);
								}
							}
							else if ((*pulse)->Type() == SpinAPI::PulseType::LongPulseStaticField)
							{

								// Create a Pulse operator in HS; only one side of exponentials as we only propagate wavevectors
								arma::sp_cx_mat pulse_operator;
								if (!space.PulseOperatorOnStatevector((*pulse), pulse_operator))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}

								// Create array containing a propagator and the current state of each system
								std::pair<arma::sp_cx_mat, arma::cx_mat> G;
								
								// Get the propagator and put it into the array together with the initial state
								arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from((A + (arma::cx_double(0.0, -1.0) * pulse_operator)) * (*pulse)->Timestep())));
								G = std::pair<arma::sp_cx_mat, arma::cx_mat>(A_sp, B);

								unsigned int steps = static_cast<unsigned int>(std::abs((*pulse)->Pulsetime() / (*pulse)->Timestep()));
								for (unsigned int n = 1; n <= steps; n++)
								{
									// Take a step, "first" is propagator and "second" is current state
									for (int i = 0; i < int(B.n_cols); i++) 
									{
										B.col(i) = G.first * G.second.col(i);
									}

									// Get the new current state vector matrix
									G.second = B;
								}
							}
							else if ((*pulse)->Type() == SpinAPI::PulseType::LongPulse)
							{
								// Create a Pulse operator in SS
								arma::sp_cx_mat pulse_operator;
								if (!space.PulseOperatorOnStatevector((*pulse), pulse_operator))
								{
									this->Log() << "Failed to create a pulse operator in SS." << std::endl;
									continue;
								}

								// Create array containing a propagator and the current state of each system
								std::pair<arma::sp_cx_mat, arma::cx_mat> G;

								// Get the propagator and put it into the array together with the initial state
								arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from((A + (arma::cx_double(0.0, -1.0) * pulse_operator * std::cos((*pulse)->Frequency() * (*pulse)->Timestep()))) * (*pulse)->Timestep())));
								G = std::pair<arma::sp_cx_mat, arma::cx_mat>(A_sp, B);

								unsigned int steps = static_cast<unsigned int>(std::abs((*pulse)->Pulsetime() / (*pulse)->Timestep()));
								for (unsigned int n = 1; n <= steps; n++)
								{
									// Take a step, "first" is propagator and "second" is current state
									for (int i = 0; i < int(B.n_cols); i++) 
									{
										B.col(i) = G.first * G.second.col(i);
									}

									// Get the new current state density vector
									G.second = B;
								}
							}
							else
							{
								this->Log() << "Not implemented yet, sorry." << std::endl;
							}

							// Get the system relax during the time

							// Create array containing a propagator and the current state of each system
							std::pair<arma::sp_cx_mat, arma::cx_mat> G;
							arma::sp_cx_mat A_sp = arma::conv_to<arma::sp_cx_mat>::from(arma::expmat(arma::conv_to<arma::cx_mat>::from(A * (*pulse)->Timestep())));
							// Get the propagator and put it into the array together with the initial state
							G = std::pair<arma::sp_cx_mat, arma::cx_mat>(A_sp, B);

							unsigned int steps = static_cast<unsigned int>(std::abs(timerelaxation / (*pulse)->Timestep()));
							for (unsigned int n = 1; n <= steps; n++)
							{
								// Take a step, "first" is propagator and "second" is current state
								for (int i = 0; i < int(B.n_cols); i++) 
								{
									B.col(i) = G.first * G.second.col(i);
								}

								// Get the new current state density vector
								G.second = B;
							}
						}
					}
				}
			}

			// Setting or calculating total time.
			double totaltime;
			this->Properties()->Get("totaltime", totaltime);

			// Setting timestep
			double dt;
			double timestep;
			this->Properties()->Get("timestep", timestep);
			dt = timestep;

			if (dt > std::pow(2, -53))
			{
				this->Log() << "Time step is chosen as " << dt << " ns." << std::endl;
			}
			else
			{
				this->Log() << "Time step is undefined. Using the default." << std::endl;
				std::cout << "# ERROR: undefined time step! Using the default of 1 ns." << std::endl;
				dt = 1;
				this->Log() << "Time step is chosen as " << dt << " ns." << std::endl;
			}

			// Number of time propagation steps
			int num_steps = std::ceil(totaltime / dt);
			this->Log() << "Number of time propagation steps: " << num_steps << "." << std::endl;
			if(num_steps == 0)
			{
				num_steps = 1;
				this->Log() << "change Number of propagation steps to: " << num_steps << " in order to propagate one step." << std::endl;
			}
			// Choose Propagation Method and other parameters
			std::string propmethod;
			this->Properties()->Get("propagationmethod", propmethod);

			std::string precision;
			this->Properties()->Get("precision", precision);

			int krylovsize;
			this->Properties()->Get("krylovsize", krylovsize);

			double krylovtol;
			this->Properties()->Get("krylovtol", krylovtol);

			if (propmethod == "autoexpm")
			{
				this->Log() << "Autoexpm is chosen as the propagation method." << std::endl;
				if (precision == "double")
				{
					this->Log() << "Double precision is chosen for the autoexpm method." << std::endl;
				}
				else if (precision == "single")
				{
					this->Log() << "Single precision is chosen for the autoexpm method." << std::endl;
				}
				else if (precision == "half")
				{
					this->Log() << "Half precision is chosen for the autoexpm method." << std::endl;
				}
				else
				{
					std::cout << "# ERROR: undefined precision. Using single digit precision!" << std::endl;
					this->Log() << "No precision for autoexpm method was defined. Using single digit precision." << std::endl;
					precision = "single";
				}
			}
			else if (propmethod == "krylov")
			{
				if (krylovsize > 0)
				{
					this->Log() << "Krylov basis size is chosen as " << krylovsize << "." << std::endl;
					if (krylovtol > 0)
					{
						this->Log() << "Tolerance for krylov propagation is chosen as " << krylovtol << "." << std::endl;
					}
					else
					{
						std::cout << "# ERROR: undefined tolerance for krylov subspace propagation! Using the default of 1e-16." << std::endl;
						this->Log() << "Undefined tolerance for the krylov subspace. Using the default of 1e-16." << std::endl;
						krylovtol = 1e-16;
					}
				}
				else
				{
					std::cout << "# ERROR: undefined size of the krylov subspace! Using the default size of 16." << std::endl;
					this->Log() << "Undefined size of the krylov subspace. Using the default size of 16." << std::endl;
					krylovsize = 16;
					if (krylovtol > 0)
					{
						this->Log() << "Tolerance for krylov propagation is chosen as " << krylovtol << "." << std::endl;
					}
					else
					{
						std::cout << "# ERROR: undefined tolerance for krylov subspace propagation! Using the default of 1e-16." << std::endl;
						this->Log() << "Undefined tolerance for the krylov subspace. Using the default of 1e-16." << std::endl;
						krylovtol = 1e-16;
					}
				}
			}
			else
			{
				std::cout << "# WARNING: Undefined propagation method, using normal exponential method."<< std::endl;  // autoexpm with single accuracy!" << std::endl;
				this->Log() << "WARNING: Undefined propagation method, using normal exponential method."<< std::endl;  // autoexpm with single accuracy." << std::endl;
				propmethod = "normal";
			}

			// Initialize time propagation placeholders
			arma::mat ExptValues;
			ExptValues.zeros(num_steps, projection_counter);
			arma::vec time(num_steps);

			// Propagate the system in time using the specified method

			// Propagation using autoexpm for matrix exponential
			if (propmethod == "autoexpm")
			{
				arma::mat M; // used for variable estimation
				// Include the recombination operator K
				H = H - arma::cx_double(0.0, 1.0) * K;
				for (int k = 0; k < num_steps; k++)
				{
					// Set the current time
					double current_time = k * dt;
					time(k) = current_time;

					// Calculate the expected values for each transition operator
					for (int idx = 0; idx < projection_counter; idx++)
					{
						double abs_trace = std::real(arma::trace(B.t() * Operators[idx] * B));
						double expected_value = abs_trace / Z;
						ExptValues(k, idx) = expected_value;
					}

					// Update B using the Higham propagator
					B = space.HighamProp(H, B, -arma::cx_double(0.0, 1.0) * dt, precision, M);
				}

				for (int k = 0; k < num_steps; k++)
				{
					// Obtain results
					this->Data() << this->RunSettings()->CurrentStep() << " ";
					this->Data() << time(k) << " ";
					this->WriteStandardOutput(this->Data());

					for (int idx = 0; idx < projection_counter; idx++)
					{
						this->Data() << " " << ExptValues(k, idx);
					}
					this->Data() << std::endl;
				}
			// Propagation using krylov subspace method for matrix exponential
			}
			else if (propmethod == "krylov")
			{
				// Include the recombination operator K
				H =  H - arma::cx_double(0.0, 1.0) * K;

				// #pragma omp parallel for
				for (int itr = 0; itr < Z; itr++)
				{
					arma::cx_vec prop_state = B.col(itr);

					// Set the current time
					double current_time = 0;
					time(0) = current_time;

					// Calculate the expected values for each transition operator
					for (int idx = 0; idx < projection_counter; idx++)
					{
						double result = std::abs(arma::cdot(prop_state, Operators[idx] * prop_state));
						ExptValues(0, idx) += result;
					}

					arma::cx_mat Hessen; // Upper Hessenberg matrix
					Hessen.zeros(krylovsize, krylovsize);

					arma::cx_mat KryBasis(InitialStateVector.n_rows * Z, krylovsize, arma::fill::zeros); // Orthogonal krylov subspace

					KryBasis.col(0) = prop_state / norm(prop_state);

					double h_mplusone_m;
					space.ArnoldiProcess(H, prop_state, KryBasis, Hessen, krylovsize, h_mplusone_m);

					arma::cx_colvec e1;
					e1.zeros(krylovsize);
					e1(0) = 1;
					arma::cx_colvec ek;
					ek.zeros(krylovsize);
					ek(krylovsize - 1) = 1;

					arma::cx_vec cx = arma::expmat(Hessen * dt) * e1;

					prop_state = norm(prop_state) * KryBasis * cx;

					int k = 1;

					while (k < num_steps)
					{
						// Set the current time
						current_time = k * dt;
						time(k) = current_time;

						// Calculate the expected values for each transition operator
						for (int idx = 0; idx < projection_counter; idx++)
						{
							double result = std::abs(arma::cdot(prop_state, Operators[idx] * prop_state));
							ExptValues(k, idx) += result;
						}

						Hessen.zeros(krylovsize, krylovsize);
						KryBasis.zeros(InitialStateVector.n_rows * Z, krylovsize);

						KryBasis.col(0) = prop_state / norm(prop_state);

						space.ArnoldiProcess(H, prop_state, KryBasis, Hessen, krylovsize, h_mplusone_m);
						cx = arma::expmat(Hessen * dt) * e1;

						// Update the state using Krylov Subspace propagator
						prop_state = norm(prop_state) * KryBasis * cx;
						k++;
					}
				}

				ExptValues /= Z;
				
				for (int k = 0; k < num_steps; k++)
				{
					// Obtain results
					this->Data() << this->RunSettings()->CurrentStep() << " ";
					this->Data() << time(k) << " ";
					this->WriteStandardOutput(this->Data());

					for (int idx = 0; idx < projection_counter; idx++)
					{
						this->Data() << " " << ExptValues(k, idx);
					}
					this->Data() << std::endl;
				}
			}
			else if (propmethod == "normal")
			{
				this->Log() << "Using robust matrix exponential propagator for time-independent Hamiltonian." << std::endl;

				// Include the recombination operator K
				arma::sp_cx_mat H_total = arma::cx_double(0.0, -1.0) * H - K;

				// Precompute the matrix exponential for the entire time step
				arma::cx_mat exp_H = arma::expmat(arma::cx_mat(H_total) * dt);

				// Propagate B
				for (int k = 0; k < num_steps; ++k) {
					// Set the current time
					double current_time = k * dt;
					time(k) = current_time;

					// Calculate the expected values for each transition operator
					for (int idx = 0; idx < projection_counter; ++idx) {
						double abs_trace = std::real(arma::trace(B.t() * arma::cx_mat(Operators[idx]) * B));
						double expected_value = abs_trace / Z;
						ExptValues(k, idx) = expected_value;
					}

					for (int i = 0; i < int(B.n_cols); i++) {
						B.col(i) = exp_H * B.col(i);
					}
				}

				// Read if the result should be integrated or not
				bool integration = false;

				if (!this->Properties()->Get("integration", integration))
				{
					this->Log() << "Failed to obtain an input for an Integtation." << std::endl;

				}
				
				if(integration==true)
				{
					// Write results
					this->Data() << this->RunSettings()->CurrentStep() << " ";
					this->WriteStandardOutput(this->Data());
					
					arma::mat ans = arma::trapz(time, ExptValues);

					for (int it = 0; it < projection_counter; it++)
					{
						// Quantym yields without correction factor
						ans(0, it) = ans(0, it) * rates(it);
						
						this->Data() << std::setprecision(6) << ans(0, it) << " ";
					}
				}
				else
				{
					for (int k = 0; k < num_steps; ++k) {
						// Write results
						this->Data() << this->RunSettings()->CurrentStep() << " ";
						this->Data() << time(k) << " ";
						this->WriteStandardOutput(this->Data());

						for (int idx = 0; idx < projection_counter; ++idx) {
							this->Data() << " " << ExptValues(k, idx);
						}
						this->Data() << std::endl;
					}
				}
			}

			this->Log() << "\nDone with SpinSystem \"" << (*i)->Name() << "\"" << std::endl;
		}
		this->Data() << std::endl;
		return true;
	}

	// Writes the header of the data file (but can also be passed to other streams)
	void TaskStaticHSDirectSpectra::WriteHeader(std::ostream &_stream)
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
	bool TaskStaticHSDirectSpectra::Validate()
	{
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

}
