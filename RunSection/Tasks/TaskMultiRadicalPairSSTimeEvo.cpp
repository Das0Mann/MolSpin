/////////////////////////////////////////////////////////////////////////
// TaskMultiRadicalPairSSTimeEvo implementation (RunSection module)

// -- Multi-system version: Allows transitions between SpinSystems --
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TaskMultiRadicalPairSSTimeEvo.h"
#include "Transition.h"
#include "Interaction.h"
#include "Settings.h"
#include "State.h"
#include "ObjectParser.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskMultiRadicalPairSSTimeEvo Constructors and Destructor
	// -----------------------------------------------------
	TaskMultiRadicalPairSSTimeEvo::TaskMultiRadicalPairSSTimeEvo(const MSDParser::ObjectParser &_parser, const RunSection &_runsection) : BasicTask(_parser, _runsection), timestep(1.0), totaltime(1.0e+4),
																																reactionOperators(SpinAPI::ReactionOperatorType::Haberkorn)
	{
	}

	TaskMultiRadicalPairSSTimeEvo::~TaskMultiRadicalPairSSTimeEvo()
	{
	}
	// -----------------------------------------------------
	// TaskMultiRadicalPairSSTimeEvo protected methods
	// -----------------------------------------------------
	bool TaskMultiRadicalPairSSTimeEvo::RunLocal()
	{
		this->Log() << "Running method StaticSS-MultiRadicalPairSystem." << std::endl;

		// If this is the first step, write first part of header to the data file
		if (this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}

		// Loop through all SpinSystems to obtain SpinSpace objects
		//auto systems = this->SpinSystems();
		auto systems = this->SpinSystems();

		if(systems.size() > 1)
		{
			this->Log() << "Error: Only supports one spin system" << std::endl;
			return false;
		}

		//get the number of sub systems
		int SubSystems = 0;
		if(!this->Properties().get()->Get("subsystems", SubSystems))
		{
			this->Log() << "Error: Number of sub systems not specified. Correct syntax is e.g. subsystems: 2" << std::endl;
			return false;
		}
		
		std::string SubSystemNames; 
		if(!this->Properties().get()->Get("subsystemnames", SubSystemNames))
		{
			this->Log() << "Note: subsystemnames not specified. Asssuming names are system1, system2, system3...." << std::endl;
		}
		
		std::vector<std::pair<std::string,std::vector<std::string>>> SubSystemSpins;
		std::vector<std::pair<std::string,std::vector<SpinAPI::interaction_ptr>>> SubSystemsInteractions;
		std::vector<std::pair<std::string,std::vector<SpinAPI::transition_ptr>>> SubSystemsTransitions;
		std::vector<std::pair<SpinAPI::transition_ptr, std::pair<std::string, std::string>>>SubSystemsInterTransitions;
		{
			std::stringstream ss(SubSystemNames);
			while(ss.good())
			{
				std::string Name;
				std::getline(ss,Name, ',');
				std::transform(Name.begin(), Name.end(), Name.begin(), [](unsigned char c){return std::tolower(c);});
				SubSystemSpins.push_back({Name, {}});
			}
		}

		//Verify that the subsytems are valid by checking if all the spin objects listed have been loaded 
		int num = 0;
		for(auto system : SubSystemSpins)
		{
			std::vector<std::string> SubSystemTemp;
			std::string SpinObjects = "";
			if(!this->Properties().get()->Get(system.first, SpinObjects))
			{
				this->Log() << "Error: Failed to Load SubSystem " << system.first << std::endl;
				std::cout << "Error: Failed to Load SubSystem " << system.first << std::endl;
				num++;
				continue;
			}

			std::stringstream ss(SpinObjects);
			while(ss.good())
			{
				std::string Name;
				std::getline(ss,Name, ',');
				SubSystemTemp.push_back(Name);
			}
			{
				auto space = std::make_shared<SpinAPI::SpinSpace>(*(*systems.cbegin()));
				bool loaded = true;
				for(std::string i : SubSystemTemp)
				{
					if(!space->Contains(i))
					{
						this->Log() << "Error: " << i << " is not a loaded spin object" << std::endl;
						std::cout << "Error: " << i << " is not a loaded spin object" << std::endl;
						loaded &= false;
						continue;
					}
					loaded &= true;
				}
				if(!loaded)
				{
					this->Log() << "Failed to load SubSystem " << system.first << std::endl;
					std::cout << "Failed to load SubSystem " << system.first << std::endl;
					num++;
					continue;
				}
				
			}

			SubSystemSpins[num].second = SubSystemTemp;
			auto FindInVec = [&SubSystemSpins](std::string str)
			{
				for(int i = 0; i < SubSystemSpins.size(); i++)
				{
					if(SubSystemSpins[i].first == str)
					{
						return true;
					}
				}
				return false;
			};

			auto MainSpinSystem = systems.cbegin();
			SubSystemsTransitions.push_back({system.first, {}});
			for(auto transiton : MainSpinSystem->get()->Transitions())
			{
				int TransitionType = 0; //0 = Transition out of one subsystem, 1 = Transition between two subsytems
				std::string SourceSubSystem = "";
				//All transitions need a source sub system
				if(transiton->Properties()->Get("sourcesubsystem", SourceSubSystem))
				{
					if(SourceSubSystem != system.first)
					{
						continue;
					}
					if(!FindInVec(SourceSubSystem))
					{
						this->Log() << "Error: No source sub system specified for transition " << transiton->Name() << std::endl;
						std::cout << "Failed to load transition " << transiton->Name() << ". No source sub system specified" << std::endl;
						transiton->SetValid(false);
						continue;
					}
				}
				else
				{
					this->Log() << "Error: No source sub system specified for transition " << transiton->Name() << std::endl;
					std::cout << "Failed to load transition " << transiton->Name() << ". No source sub system specified" << std::endl;
					transiton->SetValid(false);
					continue;
				}

				std::string TargetSubSystem = "";
				if(transiton->Properties()->Get("targetsubsystem", TargetSubSystem))
				{
					if(!FindInVec(TargetSubSystem))
					{
						this->Log() << "Error: No target sub system specified for transition " << transiton->Name() << std::endl;
						std::cout << "Failed to load transition " << transiton->Name() << ". No target sub system specified" << std::endl;
						transiton->SetValid(false);
						continue;
					}
					TransitionType = 1;
				}
				switch (TransitionType)
				{
				case 0:
					SubSystemsTransitions[num].second.push_back(transiton);
					break;
				case 1:
					SubSystemsInterTransitions.push_back({transiton, {SourceSubSystem, TargetSubSystem}});
				default:
					break;
				}
			}

			//similar system for interactions
			for(auto interaction : MainSpinSystem->get()->Interactions())
			{
				std::string SubSystem = "";
				SubSystemsInteractions.push_back({system.first, {}});
				if(!interaction->Properties()->Get("subsystem", SubSystem))
				{
					this->Log() << "Error: No sub system specified for interaction " << interaction->Name() << std::endl;
					std::cout << "Failed to load interaction " << interaction->Name() << ". No sub system specified" << std::endl;
					interaction->SetValid(false);
					continue;
				}
				if(SubSystem != system.first)
				{
					continue;
				}
				
				//TODO: checks that the spin objects involved are in that subsystem
				SubSystemsInteractions[num].second.push_back(interaction);
			}

			num++;
		}


		std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> spaces;
		unsigned int dimensions = 0;
		for (auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			auto space = std::make_shared<SpinAPI::SpinSpace>(*(*i));

			// We are using superoperator space, and need the total dimensions
			space->UseSuperoperatorSpace(true);
			space->SetReactionOperatorType(this->reactionOperators);
			dimensions += space->SpaceDimensions();

			// Make sure to save the newly created spin space
			spaces.push_back(std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>(*i, space));
		}

		// Now, create a matrix to hold the Liouvillian superoperator and the initial state
		arma::sp_cx_mat L(dimensions, dimensions);
		arma::cx_vec rho0(dimensions);
		unsigned int nextDimension = 0; // Keeps track of the dimension where the next spin space starts

		// Loop through the systems again to fill this matrix and vector
		int spinsystem = 0; 
		auto SpinSpace = spaces.cbegin();
		for(auto i = SubSystemSpins.cbegin(); i != SubSystemSpins.end(); i++)
		{
			std::cout << "Radical Pair: " << spinsystem << std::endl;
			//Get intitial state
			auto initial_states = SpinSpace->first->InitialState()[spinsystem]; //get error when the state doesn't exist, TODO: fix this
			arma::cx_mat rho0HS;
			if(initial_states->Name() == "zero")
			{
				rho0HS = arma::zeros<arma::cx_mat>(SpinSpace->second->HilbertSpaceDimensions(), SpinSpace->second->HilbertSpaceDimensions());
			}
			else
			{
				if(!SpinSpace->second->GetState(initial_states, rho0HS))
				{
					this->Log() << "ERROR: Failed to obtain projection matrix onto state \"" << initial_states->Name() << "\", initial state of SpinSubSystem \"" << i->first << "\"." << std::endl;
					return false;
				}
				rho0HS /= arma::trace(rho0HS);
			}
			spinsystem++;
		}
		for (auto i = spaces.cbegin(); i != spaces.cend(); i++)
		{
			std::cout << spinsystem << std::endl;
			// Make sure we have an initial state
			auto initial_states = i->first->InitialState();
			arma::cx_mat rho0HS;
			if (initial_states.size() < 1)
			{
				this->Log() << "Note: No initial state specified for spin system \"" << i->first->Name() << "\", setting the initial state to zero." << std::endl;
				rho0HS = arma::zeros<arma::cx_mat>(i->second->HilbertSpaceDimensions(), i->second->HilbertSpaceDimensions());
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
						rho0HS = tmp_rho0;
					else
						rho0HS += tmp_rho0;
				}
				rho0HS /= arma::trace(rho0HS); // The density operator should have a trace of 1
			}
			// Now put the initial state into the superspace vector
			arma::cx_vec rho0vec;
			if (!i->second->OperatorToSuperspace(rho0HS, rho0vec))
			{
				this->Log() << "ERROR: Failed convert initial state to superspace for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}
			rho0.rows(nextDimension, nextDimension + i->second->SpaceDimensions() - 1) = rho0vec;

			// Next, get the Hamiltonian
			arma::sp_cx_mat H;
			if (!i->second->Hamiltonian(H))
			{
				this->Log() << "ERROR: Failed to obtain the superspace Hamiltonian for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}
			L.submat(nextDimension, nextDimension, nextDimension + i->second->SpaceDimensions() - 1, nextDimension + i->second->SpaceDimensions() - 1) = arma::cx_double(0.0, -1.0) * H;

			// Then get the reaction operators
			arma::sp_cx_mat K;
			if (!i->second->TotalReactionOperator(K))
			{
				this->Log() << "ERROR: Failed to obtain matrix representation of the reaction operators for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}

			L.submat(nextDimension, nextDimension, nextDimension + i->second->SpaceDimensions() - 1, nextDimension + i->second->SpaceDimensions() - 1) -= K;

			//this entire system needs to go
			// Obtain the creation operators - note that we need to loop through the other SpinSystems again to find transitions leading into the current SpinSystem
			unsigned int nextCDimension = 0; // Similar to nextDimension, but to keep track of first dimension for this other SpinSystem
			for (auto j = spaces.cbegin(); j != spaces.cend(); j++)
			{
				// Creation operators are off-diagonal elements
				if (j != i)
				{
					// Check all transitions whether they should produce a creation operator
					for (auto t = j->first->transitions_cbegin(); t != j->first->transitions_cend(); t++)
					{
						// Does the Transition lead into the current spin space?
						if ((*t)->Target() == i->first)
						{
							// Prepare a creation operator
							arma::sp_cx_mat C;
							if(!j->second->ReactionOperator((*t), C))
							{
								this->Log() << "ERROR: Failed to obtain matrix representation of the  reaction operator for transition \"" << (*t)->Name() << "\"!" << std::endl;
								return false;
							}

							// Put it into the total Liouvillian:
							//  - The row should be that of the current spin space (the target space)
							//  - The column should be that of the source spin space (the spin system containing the Transition object)
							L.submat(nextDimension, nextCDimension, nextDimension + i->second->SpaceDimensions() - 1, nextCDimension + j->second->SpaceDimensions() - 1) += C;
						}
					}
				}

				// Move on to check next spin system for transitions into the current spin space
				nextCDimension += j->second->SpaceDimensions();
			}

			// Move on to next spin space
			nextDimension += i->second->SpaceDimensions();
			spinsystem++;
		}

		// Write results for initial state as well (i.e. at time 0)
		this->Data() << this->RunSettings()->CurrentStep() << " 0 ";
		this->WriteStandardOutput(this->Data());
		nextDimension = 0;
		for (auto i = spaces.cbegin(); i != spaces.cend(); i++)
		{
			// Get the superspace result vector and convert it back to the native Hilbert space
			arma::cx_mat rho_result;
			arma::cx_vec rho_result_vec;
			rho_result_vec = rho0.rows(nextDimension, nextDimension + i->second->SpaceDimensions() - 1);
			if (!i->second->OperatorFromSuperspace(rho_result_vec, rho_result))
			{
				this->Log() << "ERROR: Failed to convert resulting superspace-vector back to native Hilbert space for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}

			// Get the results
			this->GatherResults(rho_result, *(i->first), *(i->second));

			// Move on to next spin space
			nextDimension += i->second->SpaceDimensions();
		}
		this->Data() << std::endl;

		// We need the propagator
		this->Log() << "Calculating the propagator..." << std::endl;
		arma::cx_mat P = arma::expmat(arma::conv_to<arma::cx_mat>::from(L) * this->timestep);

		// Perform the calculation
		this->Log() << "Ready to perform calculation." << std::endl;
		unsigned int steps = static_cast<unsigned int>(std::abs(this->totaltime / this->timestep));
		for (unsigned int n = 1; n <= steps; n++)
		{
			// Write first part of the data output
			this->Data() << this->RunSettings()->CurrentStep() << " ";
			this->Data() << (static_cast<double>(n) * this->timestep) << " ";
			this->WriteStandardOutput(this->Data());

			// Propagate (use special scope to be able to dispose of the temporary vector asap)
			{
				arma::cx_vec tmp = P * rho0;
				rho0 = tmp;
			}

			// Retrieve the resulting density matrix for each spin system and output the results
			nextDimension = 0;
			for (auto i = spaces.cbegin(); i != spaces.cend(); i++)
			{
				// Get the superspace result vector and convert it back to the native Hilbert space
				arma::cx_mat rho_result;
				arma::cx_vec rho_result_vec;
				rho_result_vec = rho0.rows(nextDimension, nextDimension + i->second->SpaceDimensions() - 1);
				if (!i->second->OperatorFromSuperspace(rho_result_vec, rho_result))
				{
					this->Log() << "ERROR: Failed to convert resulting superspace-vector back to native Hilbert space for spin system \"" << i->first->Name() << "\"!" << std::endl;
					return false;
				}

				// Get the results
				this->GatherResults(rho_result, *(i->first), *(i->second));

				// Move on to next spin space
				nextDimension += i->second->SpaceDimensions();
			}

			// Terminate the line in the data file after iteration through all spin systems
			this->Data() << std::endl;
		}

		this->Log() << "Done with calculation." << std::endl;

		return true;
	}

	// Gathers and outputs the results from a given time-integrated density operator
	void TaskMultiRadicalPairSSTimeEvo::GatherResults(const arma::cx_mat &_rho, const SpinAPI::SpinSystem &_system, const SpinAPI::SpinSpace &_space)
	{
		// Loop through all states
		arma::cx_mat P;
		auto states = _system.States();
		for (auto j = states.cbegin(); j != states.cend(); j++)
		{
			if (!_space.GetState((*j), P))
			{
				this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << _system.Name() << "\"." << std::endl;
				continue;
			}

			// Return the yield for this state - note that no reaction rates are included here.
			this->Data() << std::abs(arma::trace(P * _rho)) << " ";
		}
	}

	// Writes the header of the data file (but can also be passed to other streams)
	void TaskMultiRadicalPairSSTimeEvo::WriteHeader(std::ostream &_stream)
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

	// Validation of the required input
	bool TaskMultiRadicalPairSSTimeEvo::Validate()
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
	// -----------------------------------------------------
}
