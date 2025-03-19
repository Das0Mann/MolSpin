/////////////////////////////////////////////////////////////////////////
// TaskMultiDynamicHSTimeEvo implementation (RunSection module)
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TaskMultiDynamicHSTimeEvo.h"
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "SpinSystem.h"
#include "Interaction.h"
#include "ObjectParser.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskStaticSS Constructors and Destructor
	// -----------------------------------------------------
	TaskMultiDynamicHSTimeEvo::TaskMultiDynamicHSTimeEvo(const MSDParser::ObjectParser &_parser, const RunSection &_runsection) : BasicTask(_parser, _runsection),
																																  timestep(0.01), totaltime(1.0e+4), outputstride(1)
	{
	}

	TaskMultiDynamicHSTimeEvo::~TaskMultiDynamicHSTimeEvo()
	{
	}
	// -----------------------------------------------------
	// TaskStaticSS protected methods
	// -----------------------------------------------------
	bool TaskMultiDynamicHSTimeEvo::RunLocal()
	{
		this->Log() << "Running method DynamicHS-MultiSystem." << std::endl;

		// If this is the first step, write first part of header to the data file
		if (this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}

		// Loop through all SpinSystems to obtain SpinSpace objects
		auto systems = this->SpinSystems();
		std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> spaces;
		for (auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Create a SpinSpace for the spin system
			auto space = std::make_shared<SpinAPI::SpinSpace>(*(*i));
			space->UseSuperoperatorSpace(false);
			space->SetTime(0.0);
			spaces.push_back(std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>(*i, space));
		}

		// We need a Hamiltonian, total reaction operator, creation operators, and a current state for each spin system
		std::vector<arma::cx_mat> H;
		std::vector<arma::cx_mat> K;
		std::vector<arma::cx_mat> rho;
		std::vector<arma::sp_cx_mat> C; // Reaction operators in target space (produced by creation operators acting on density operators, or similar means). Note: Sparse matrices!

		// Dynamic parts
		std::vector<arma::sp_cx_mat> dH;
		std::vector<arma::sp_cx_mat> dK;

		bool timedependentInteractions = false;
		bool timedependentTransitions = false;

		// Loop through the systems again to fill these matrices
		for (auto i = spaces.cbegin(); i != spaces.cend(); i++)
		{
			// Make sure we have an initial state
			auto initial_states = i->first->InitialState();
			arma::cx_mat rho0;
			if (initial_states.size() < 1)
			{
				this->Log() << "Note: No initial state specified for spin system \"" << i->first->Name() << "\", setting the initial state to zero." << std::endl;
				rho0 = arma::zeros<arma::cx_mat>(i->second->SpaceDimensions(), i->second->SpaceDimensions());
			}
			else
			{
				// Get the initial state for the current system
				for (auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
				{
					arma::cx_mat tmp_rho0;
					if (!i->second->GetState(*j, tmp_rho0))
					{
						this->Log() << "ERROR: Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << i->first->Name() << "\"." << std::endl;
						return false;
					}
					if (j == initial_states.cbegin())
						rho0 = tmp_rho0;
					else
						rho0 += tmp_rho0;
				}
				rho0 /= arma::trace(rho0); // The density operator should have a trace of 1
			}

			// Put the created initial state density operator into the list
			rho.push_back(rho0);

			// Next, get the Hamiltonian
			arma::cx_mat tmp_H;
			arma::sp_cx_mat tmp_dH = arma::sp_cx_mat(i->second->SpaceDimensions(), i->second->SpaceDimensions());
			if (!i->second->StaticHamiltonian(tmp_H) || (i->second->HasTimedependentInteractions() && !i->second->DynamicHamiltonian(tmp_dH)))
			{
				this->Log() << "ERROR: Failed to obtain the Hamiltonian for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}
			H.push_back(arma::cx_double(0.0, -1.0) * tmp_H);
			dH.push_back(arma::cx_double(0.0, -1.0) * tmp_dH);

			// Then get the reaction operators (decay part)
			arma::cx_mat tmp_K;
			arma::sp_cx_mat tmp_dK = arma::sp_cx_mat(i->second->SpaceDimensions(), i->second->SpaceDimensions());
			if (!i->second->StaticTotalReactionOperator(tmp_K) || (i->second->HasTimedependentTransitions() && !i->second->DynamicTotalReactionOperator(tmp_dK)))
			{
				this->Log() << "ERROR: Failed to obtain matrix representation of the reaction operators for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}
			K.push_back(tmp_K);
			dK.push_back(tmp_dK);

			// Allocate matrices for the creation operators
			C.push_back(arma::sp_cx_mat(i->second->SpaceDimensions(), i->second->SpaceDimensions()));

			// Check whether we have any time-dependent interactions or transitions
			timedependentInteractions |= i->second->HasTimedependentInteractions();
			timedependentTransitions |= i->second->HasTimedependentTransitions();
		}

		// ---------------------------------------------------------------------------------------------------
		// Perform the calculation
		// ---------------------------------------------------------------------------------------------------
		this->Log() << "Ready to perform calculation." << std::endl;
		this->OutputResults(spaces, rho, 0); // Write initial state results
		unsigned int steps = static_cast<unsigned int>(std::abs(this->totaltime / this->timestep));
		for (unsigned int n = 0; n < steps; n++)
		{
			// Obtain creation operators
			this->GetCreationOperators(spaces, C, rho);

			// Propagate
			// this->AdvanceStep_AsyncLeapfrog(spaces, H, dH, K, dK, C, rho);
			this->AdvanceStep_RK4(spaces, H, dH, K, dK, C, rho);

			// Print results
			if (n % this->outputstride == 0)
				this->OutputResults(spaces, rho, n + 1);

			// Update time-dependent interactions or reaction operators
			if ((timedependentInteractions || timedependentTransitions) && n < steps - 1)
				this->UpdateTimeDependences(spaces, dH, dK, n + 1);
		}
		// ---------------------------------------------------------------------------------------------------

		this->Log() << "Done with calculation." << std::endl;

		return true;
	}

	// -----------------------------------------------------
	// The method that prepares all the creation operators
	void TaskMultiDynamicHSTimeEvo::GetCreationOperators(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> &_spaces,
														 std::vector<arma::sp_cx_mat> &_C, const std::vector<arma::cx_mat> &_rho)
	{
		// Reset all creation operators
		for (unsigned int i = 0; i < _spaces.size(); i++)
			_C[i] = arma::sp_cx_mat(_spaces[i].second->SpaceDimensions(), _spaces[i].second->SpaceDimensions());

		// Obtain the new creation operators
		for (unsigned int i = 0; i < _spaces.size(); i++)
		{
			// Loop through all transitions in each spin system
			for (auto j = _spaces[i].first->transitions_cbegin(); j != _spaces[i].first->transitions_cend(); j++)
			{
				// If the transition has a target spin system
				if ((*j)->Target() != nullptr && (*j)->TargetState() != nullptr)
				{
					// Find the index of the target system
					auto target = _spaces.size();
					for (unsigned int n = 0; n < _spaces.size(); n++)
					{
						if (_spaces[n].first == (*j)->Target())
						{
							target = n;
							break;
						}
					}

					// Check that the target system was found
					if (target >= _spaces.size())
					{
						this->Log() << "ERROR: Failed to properly generate creation operators." << std::endl;
						return;
					}

					// Get the creation rate
					arma::cx_mat P;
					if (!_spaces[i].second->GetState((*j)->SourceState(), P))
					{
						this->Log() << "Failed to obtain projection matrix onto source state of transition \"" << (*j)->Name() << "\"." << std::endl;
						continue;
					}
					double cRate = std::abs(arma::trace(P * _rho[i]));

					// Generate the creation operator
					arma::sp_cx_mat sP;
					if (!_spaces[target].second->ReactionTargetOperator((*j), cRate, sP))
					{
						this->Log() << "Failed to obtain reaction operator for target state of transition \"" << (*j)->Name() << "\"." << std::endl;
						continue;
					}
					_C[target] += sP;
				}
			}
		}
	}

	// -----------------------------------------------------
	// The timestep function
	void TaskMultiDynamicHSTimeEvo::AdvanceStep_AsyncLeapfrog(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> &_spaces,
															  const std::vector<arma::cx_mat> &_H, const std::vector<arma::sp_cx_mat> &_dH,
															  const std::vector<arma::cx_mat> &_K, const std::vector<arma::sp_cx_mat> &_dK,
															  const std::vector<arma::sp_cx_mat> &_C, std::vector<arma::cx_mat> &_rho)
	{
		// Loop through all spaces and propagate each individually
		for (unsigned int i = 0; i < _spaces.size(); i++)
		{
			arma::cx_mat dRho = arma::cx_double(0.0, -1.0) * ((_H[i] + _dH[i]) * _rho[i] - _rho[i] * (_H[i] + _dH[i])) - ((_K[i] + _dK[i]) * _rho[i] + _rho[i] * (_K[i] + _dK[i])) + _C[i];
			_rho[i] += dRho * this->timestep / 2.0;
			dRho = (arma::cx_double(0.0, -1.0) * ((_H[i] + _dH[i]) * _rho[i] - _rho[i] * (_H[i] + _dH[i])) - ((_K[i] + _dK[i]) * _rho[i] + _rho[i] * (_K[i] + _dK[i])) + _C[i]);
			_rho[i] += dRho * this->timestep / 2.0;
		}
	}

	void TaskMultiDynamicHSTimeEvo::AdvanceStep_RK4(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> &_spaces,
													const std::vector<arma::cx_mat> &_H, const std::vector<arma::sp_cx_mat> &_dH,
													const std::vector<arma::cx_mat> &_K, const std::vector<arma::sp_cx_mat> &_dK,
													const std::vector<arma::sp_cx_mat> &_C, std::vector<arma::cx_mat> &_rho)
	{
		for (unsigned int i = 0; i < _spaces.size(); i++)
		{
			arma::cx_mat k1 = ComputeRhoDot(_H[i], _dH[i], _K[i], _dK[i], _C[i], _rho[i]);
			arma::cx_mat k2 = ComputeRhoDot(_H[i], _dH[i], _K[i], _dK[i], _C[i], _rho[i] + 0.5 * this->timestep * k1);
			arma::cx_mat k3 = ComputeRhoDot(_H[i], _dH[i], _K[i], _dK[i], _C[i], _rho[i] + 0.5 * this->timestep * k2);
			arma::cx_mat k4 = ComputeRhoDot(_H[i], _dH[i], _K[i], _dK[i], _C[i], _rho[i] + this->timestep * k3);

			_rho[i] += (this->timestep / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
		}
	}

	arma::cx_mat TaskMultiDynamicHSTimeEvo::ComputeRhoDot(const arma::cx_mat &H, const arma::sp_cx_mat &dH,
														  const arma::cx_mat &K, const arma::sp_cx_mat &dK,
														  const arma::sp_cx_mat &C, const arma::cx_mat &rho)
	{
		return arma::cx_double(0.0, -1.0) * ((H + dH) * rho - rho * (H + dH)) - ((K + dK) * rho + rho * (K + dK)) + C;
	}

	// -----------------------------------------------------
	// Writes the output for a timestep
	void TaskMultiDynamicHSTimeEvo::OutputResults(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> &_spaces,
												  const std::vector<arma::cx_mat> &_rho, const unsigned int _step)
	{
		// Obtain the results
		arma::cx_mat PState;
		this->Data() << this->RunSettings()->CurrentStep() << " ";
		this->Data() << (static_cast<double>(_step) * this->timestep) << " ";
		this->WriteStandardOutput(this->Data());
		for (unsigned int i = 0; i < _spaces.size(); i++)
		{
			auto states = _spaces[i].first->States();
			for (auto j = states.cbegin(); j != states.cend(); j++)
			{
				if (!_spaces[i].second->GetState((*j), PState))
				{
					this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\"." << std::endl;
					continue;
				}

				this->Data() << std::abs(arma::trace(PState * _rho[i])) << " ";
			}
		}

		// Terminate the line in the data file after iteration through all steps
		this->Data() << std::endl;
	}

	// -----------------------------------------------------
	// Method to update time-dependent interactions and reactions
	void TaskMultiDynamicHSTimeEvo::UpdateTimeDependences(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> &_spaces,
														  std::vector<arma::sp_cx_mat> &_H, std::vector<arma::sp_cx_mat> &_K, unsigned int _step)
	{
		// Handle each spin space
		for (unsigned int i = 0; i < _spaces.size(); i++)
		{
			// First update the time
			_spaces[i].second->SetTime(static_cast<double>(_step) * this->timestep);

			// Do we need to update the interactions?
			if (_spaces[i].second->HasTimedependentInteractions())
			{
				if (!_spaces[i].second->DynamicHamiltonian(_H[i]))
				{
					this->Log() << "ERROR: Failed to update the Hamiltonian for spin system \"" << _spaces[i].first->Name() << "\"!" << std::endl;
					return;
				}
				_H[i] *= arma::cx_double(0.0, -1.0);
			}

			// Do we need to update the transitions?
			if (_spaces[i].second->HasTimedependentTransitions())
			{
				if (!_spaces[i].second->DynamicTotalReactionOperator(_K[i]))
				{
					this->Log() << "ERROR: Failed to update matrix representation of the reaction operators for spin system \"" << _spaces[i].first->Name() << "\"!" << std::endl;
					return;
				}
			}
		}
	}

	// -----------------------------------------------------
	// Writes the header of the data file (but can also be passed to other streams)
	void TaskMultiDynamicHSTimeEvo::WriteHeader(std::ostream &_stream)
	{
		_stream << "Step ";
		_stream << "Time(ns) ";
		this->WriteStandardOutputHeader(_stream);

		// Get header for each spin system
		auto systems = this->SpinSystems();
		for (auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Write each state name
			auto states = (*i)->States();
			for (auto j = states.cbegin(); j != states.cend(); j++)
				_stream << (*i)->Name() << "." << (*j)->Name() << " ";
		}
		_stream << std::endl;
	}

	// -----------------------------------------------------
	// Task validation method
	bool TaskMultiDynamicHSTimeEvo::Validate()
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
				// We can run the calculation if an invalid timestep was specified
				return false;
			}
		}

		// Get totaltime
		if (this->Properties()->Get("totaltime", inputTotaltime))
		{
			if (std::isfinite(inputTotaltime) && inputTotaltime > 0.0)
			{
				this->totaltime = inputTotaltime;
			}
			else
			{
				// We can run the calculation if an invalid total time was specified
				return false;
			}
		}

		// Get the output stride (i.e. write output for any n'th step, where n is the stride) - Only available for TimeEvolution calculations!
		if (this->Properties()->Get("outputstride", this->outputstride))
		{
			if (this->outputstride == 0)
			{
				this->Log() << "Cannot use outputstride 0, must be a non-zero positive integer! Using default value of 1." << std::endl;
				this->outputstride = 1;
			}
		}

		return true;
	}
	// -----------------------------------------------------
}
