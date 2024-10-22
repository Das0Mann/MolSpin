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
#include "SubSystem.h"


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
		auto systems = this->SpinSystems();
		if(systems.size() > 1)
		{
			this->Log() << "Error: Only supports one spin system" << std::endl;
			return false;
		}

		//get the number of sub systems
		auto ss = systems[0]->SubSystems();
		int SubSystems = ss.size();

		//change to this->Properties()->Get()
		//if(!this->Properties().get()->Get("subsystems", SubSystems))
		//{
		//	this->Log() << "Error: Number of sub systems not specified. Correct syntax is e.g. subsystems: 2" << std::endl;
		//	return false;
		//}
		
		//std::string SubSystemNames; 
		//if(!this->Properties().get()->Get("subsystemnames", SubSystemNames))
		//{
		//	this->Log() << "Note: subsystemnames not specified. Asssuming names are system1, system2, system3...." << std::endl;
		//}
		
		std::vector<std::pair<std::string,std::vector<std::string>>> SubSystemSpins;
		std::vector<std::pair<std::string,std::vector<SpinAPI::interaction_ptr>>> SubSystemsInteractions;
		std::vector<SubSystemTransition> SubSystemsTransitions;
		//std::vector<std::pair<SpinAPI::transition_ptr, std::pair<std::string, std::string>>>SubSystemsInterTransitions;
		//{
		//	std::stringstream ss(SubSystemNames);
		//	while(ss.good())
		//	{
		//		std::string Name;
		//		std::getline(ss,Name, ',');
		//		std::transform(Name.begin(), Name.end(), Name.begin(), [](unsigned char c){return std::tolower(c);});
		//		SubSystemSpins.push_back({Name, {}});
		//	}
		//}

		int num = 0;
		for(auto s : ss)
		{
			SubSystemSpins.push_back({s->Name(), s->GetSpinNames()});

			for(auto tr : s->GetTransitions())
			{
				//remove duplicates
				bool duplicate = false;
				for(int i = 0; i < SubSystemsTransitions.size(); i++)
				{
					if(*tr.TransitionObject == SubSystemsTransitions[i].transition)
					{
						std::string d = "";
						SubSystemsTransitions[i].transition->Properties()->Get("duplicate", d);
						if(d[0] != 't') //only need to check first letter
							duplicate = true;
							break;
					}

				}

				if(duplicate) {continue;}

				SubSystemTransition TransitionData;
				TransitionData.transition = *tr.TransitionObject;
				TransitionData.type = tr.type;
				TransitionData.source = tr.source;

				switch(tr.type)
				{
				case 0:
					SubSystemsTransitions.push_back(TransitionData);
					break;
				case 1:
					TransitionData.target = tr.target;
					SubSystemsTransitions.push_back(TransitionData);
					break;
				default:
					break;
				}
			}

			SubSystemsInteractions.push_back({s->Name(), s->GetInteractions()});

		}

		//Verify that the subsytems are valid by checking if all the spin objects listed have been loaded 
		//num = 0;
		//for(auto system : SubSystemSpins)
		//{
		//	//std::vector<std::string> SubSystemTemp;
		//	//std::string SpinObjects = "";
		//	//if(!this->Properties().get()->Get(system.first, SpinObjects))
		//	//{
		//	//	this->Log() << "Error: Failed to Load SubSystem " << system.first << std::endl;
		//	//	std::cout << "Error: Failed to Load SubSystem " << system.first << std::endl;
		//	//	num++;
		//	//	continue;
		//	//}

		//	//std::stringstream ss(SpinObjects);
		//	//while(ss.good())
		//	//{
		//	//	std::string Name;
		//	//	std::getline(ss,Name, ',');
		//	//	SubSystemTemp.push_back(Name);
		//	//}
		//	//{
		//	//	auto space = std::make_shared<SpinAPI::SpinSpace>(*(*systems.cbegin()));
		//	//	bool loaded = true;
		//	//	for(std::string i : SubSystemTemp)
		//	//	{
		//	//		if(!space->Contains(i))
		//	//		{
		//	//			this->Log() << "Error: " << i << " is not a loaded spin object" << std::endl;
		//	//			std::cout << "Error: " << i << " is not a loaded spin object" << std::endl;
		//	//			loaded &= false;
		//	//			continue;
		//	//		}
		//	//		loaded &= true;
		//	//	}
		//	//	if(!loaded)
		//	//	{
		//	//		this->Log() << "Failed to load SubSystem " << system.first << std::endl;
		//	//		std::cout << "Failed to load SubSystem " << system.first << std::endl;
		//	//		num++;
		//	//		continue;
		//	//	}
		//	//	
		//	//}

		//	//SubSystemSpins[num].second = SubSystemTemp;
		//	auto FindInVec = [&SubSystemSpins](std::string str)
		//	{
		//		for(int i = 0; i < SubSystemSpins.size(); i++)
		//		{
		//			if(SubSystemSpins[i].first == str)
		//			{
		//				return true;
		//			}
		//		}
		//		return false;
		//	};

		//	auto MainSpinSystem = systems.cbegin();

		//	//SubSystemsTransitions.push_back({system.first, {}});
		//	for(auto transiton : MainSpinSystem->get()->Transitions())
		//	{
		//		int TransitionType = 0; //0 = Transition out of one subsystem, 1 = Transition between two subsytems
		//		std::string SourceSubSystem = "";
		//		//All transitions need a source sub system
		//		if(transiton->Properties()->Get("sourcesubsystem", SourceSubSystem))
		//		{
		//			if(SourceSubSystem != system.first)
		//			{
		//				continue;
		//			}
		//			if(!FindInVec(SourceSubSystem))
		//			{
		//				this->Log() << "Error: No source sub system specified for transition " << transiton->Name() << std::endl;
		//				std::cout << "Failed to load transition " << transiton->Name() << ". No source sub system specified" << std::endl;
		//				transiton->SetValid(false);
		//				continue;
		//			}
		//		}
		//		else
		//		{
		//			this->Log() << "Error: No source sub system specified for transition " << transiton->Name() << std::endl;
		//			std::cout << "Failed to load transition " << transiton->Name() << ". No source sub system specified" << std::endl;
		//			transiton->SetValid(false);
		//			continue;
		//		}

		//		std::string TargetSubSystem = "";
		//		if(transiton->Properties()->Get("targetsubsystem", TargetSubSystem))
		//		{
		//			if(!FindInVec(TargetSubSystem))
		//			{
		//				this->Log() << "Error: No target sub system specified for transition " << transiton->Name() << std::endl;
		//				std::cout << "Failed to load transition " << transiton->Name() << ". No target sub system specified" << std::endl;
		//				transiton->SetValid(false);
		//				continue;
		//			}
		//			TransitionType = 1;
		//		}

		//		SubSystemTransition TransitionData;
		//		TransitionData.transition = transiton;
		//		TransitionData.type = TransitionType;
		//		TransitionData.source = SourceSubSystem;
		//		TransitionData.type = TransitionType;

		//		switch (TransitionType)
		//		{
		//		case 0:
		//			SubSystemsTransitions.push_back(TransitionData);
		//			break;
		//		case 1:
		//			TransitionData.target = TargetSubSystem;
		//			SubSystemsTransitions.push_back(TransitionData);
		//			break;
		//		default:
		//			break;
		//		}
		//	}

		//	//similar system for interactions
		//	for(auto interaction : MainSpinSystem->get()->Interactions())
		//	{
		//		std::string SubSystem = "";
		//		SubSystemsInteractions.push_back({system.first, {}});
		//		if(!interaction->Properties()->Get("subsystem", SubSystem))
		//		{
		//			this->Log() << "Error: No sub system specified for interaction " << interaction->Name() << std::endl;
		//			std::cout << "Failed to load interaction " << interaction->Name() << ". No sub system specified" << std::endl;
		//			interaction->SetValid(false);
		//			continue;
		//		}

		//		std::vector<std::string> systems; 
		//		std::stringstream ss(SubSystem);
		//		while(ss.good())
		//		{
		//			std::string Name;
		//			std::getline(ss,Name, ',');
		//			systems.push_back(Name);
		//		}
		//		
		//		for(auto i : systems)
		//		{
		//			if(i != system.first)
		//			{
		//				continue;
		//			}

		//			//TODO: checks that the spin objects involved are in that subsystem
		//			SubSystemsInteractions[num].second.push_back(interaction);
		//			break;
		//		}
		//	}

		//	num++;
		//}

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
		int super_dimensions = dimensions * SubSystems; //every SubSystem has the same dimension so the total operator dimension has to be bigger depending on the number of subsystems
		// Now, create a matrix to hold the Liouvillian superoperator and the initial state
		arma::sp_cx_mat L(super_dimensions, super_dimensions);
		arma::cx_vec rho0(super_dimensions);
		unsigned int nextDimension = 0; // Keeps track of the dimension where the next spin space starts

		// Loop through the systems again to fill this matrix and vector
		int spinsystem = 0; 
		auto SpinSpace = spaces.cbegin();

		std::vector<SpinAPI::state_ptr> AllStates;
		std::vector<std::string> InitialStates;
		AllStates = SpinSpace->first->GetAllStates();
		std::string InitialStatesString;
		if(!SpinSpace->first->GetProperties()->Get("initialstate", InitialStatesString))
		{
			this->Log() << "No initial states provided, density matrix will be set to 0" << std::endl;
			for(int i = 0; i < SubSystems; i++)
			{
				InitialStates.push_back({""});
			}
		}
		else
		{
			std::stringstream ss(InitialStatesString);
			while(ss.good())
			{
				std::string state;
				std::getline(ss,state, ',');
				InitialStates.push_back(state);
			}
		}
		
		for(auto i = SubSystemSpins.cbegin(); i != SubSystemSpins.end(); i++)
		{
			std::cout << "Radical Pair: " << spinsystem << std::endl;
			//Get intitial state
			auto initial_state = InitialStates[spinsystem]; 
			auto FindState = [&initial_state](auto state)
			{
				return state->Name() == initial_state;
			};
			bool found = true;
			auto FoundState = std::find_if(AllStates.begin(), AllStates.end(), FindState);
			if(*FoundState == nullptr)
			{
				found = false;
			}
			arma::cx_mat rho0HS;
			if(initial_state == "zero" || initial_state == "")
			{
				rho0HS = arma::zeros<arma::cx_mat>(SpinSpace->second->HilbertSpaceDimensions(), SpinSpace->second->HilbertSpaceDimensions());
			}
			else if(found == true)
			{
				if(!SpinSpace->second->GetState((*FoundState), rho0HS))
				{
					this->Log() << "ERROR: Failed to obtain projection matrix onto state \"" << (*FoundState)->Name() << "\", initial state of SpinSubSystem \"" << i->first << "\"." << std::endl;
					return false;
				}
				rho0HS /= arma::trace(rho0HS);
			}
			else
			{
				rho0HS = arma::zeros<arma::cx_mat>(SpinSpace->second->HilbertSpaceDimensions(), SpinSpace->second->HilbertSpaceDimensions());
			}

			// Now put the initial state into the superspace vector
			arma::cx_vec rho0vec;
			if (!SpinSpace->second->OperatorToSuperspace(rho0HS, rho0vec))
			{
				this->Log() << "ERROR: Failed convert initial state to superspace for spin subsystem \"" << SubSystemSpins[spinsystem].first << "\"!" << std::endl;
				return false;
			}
			rho0.rows(nextDimension, nextDimension + SpinSpace->second->SpaceDimensions() - 1) = rho0vec;

			// Next, get the Hamiltonian
			arma::sp_cx_mat H;
			if (!GenerateHamiltonian(SubSystemsInteractions[spinsystem].second, H, dimensions, SpinSpace->second))
			{
				this->Log() << "ERROR: Failed to obtain the superspace Hamiltonian for spin spin subsystem \"" << SubSystemSpins[spinsystem].first << "\"!" << std::endl;
				return false;
			}
			auto offset = SpinSpace->second->SpaceDimensions();
			L.submat(nextDimension, nextDimension, nextDimension + SpinSpace->second->SpaceDimensions() - 1, nextDimension + SpinSpace->second->SpaceDimensions() - 1) = arma::cx_double(0.0, -1.0) * H;
			
			// Then get the reaction operators
			arma::sp_cx_mat K;
			if (!GenerateReactionOperator(SubSystemsTransitions, K, dimensions, SpinSpace->second, i->first))
			{
				this->Log() << "ERROR: Failed to obtain matrix representation of the reaction operators for spin spin subsystem \"" << SubSystemSpins[spinsystem].first << "\"!" << std::endl;
				return false;
			}
			L.submat(nextDimension, nextDimension, nextDimension + SpinSpace->second->SpaceDimensions() - 1, nextDimension + SpinSpace->second->SpaceDimensions() - 1) -= K;

			//Loop through all the interactions that involve transitions to and from other subsystems
			unsigned int nextCDimension = 0; // Similar to nextDimension, but to keep track of first dimension for this other SpinSystem
			for (auto j = SubSystemsTransitions.begin(); j != SubSystemsTransitions.end(); j++)
			{
				//only process the transition if the transition is leading into the current spin state
				if(j->type == 0)
				{
					continue;
				}
				if(j->target != SubSystemSpins[spinsystem].first)
				{
					continue;
				}

				//identify what column the source subsystem is
				std::string SubSystemName = j->source;
				auto FindSubSystem = [&SubSystemName](std::pair<std::string,std::vector<std::string>> t)
				{
					return SubSystemName == t.first;		
				};
				auto temp = std::find_if(SubSystemSpins.begin(), SubSystemSpins.end(), FindSubSystem);
				int col = temp - SubSystemSpins.begin();
				nextCDimension = col * dimensions;

				//transition is leading into subsystem
				arma::sp_cx_mat T;
				if(!SpinSpace->second->ReactionOperator(j->transition, T))
				{
					this->Log() << "ERROR: Failed to obtain matrix representation of the reaction operator for transition \"" << j->transition->Name() << "\"!" << std::endl;
					return false;
				}
				//std::cout << T << std::endl;
				L.submat(nextDimension, nextCDimension, nextDimension + SpinSpace->second->SpaceDimensions() - 1, nextCDimension + SpinSpace->second->SpaceDimensions() - 1) += T;
				//T.print();
				//modify L for the transition leading out of the subsytem;
				//L.submat(nextCDimension, nextCDimension, nextCDimension + SpinSpace->second->SpaceDimensions() - 1, nextCDimension + SpinSpace->second->SpaceDimensions() - 1) += -1 * T;
				//std::cout << L << std::endl;
			}
			spinsystem++;
			nextDimension += dimensions;
		}
		// Write results for initial state as well (i.e. at time 0)
		this->Data() << this->RunSettings()->CurrentStep() << " 0 ";
		this->WriteStandardOutput(this->Data());
		nextDimension = 0;

		//modify this so it's the same that I've got in my code, using the vectors and superspace rather than matracies;

		struct TrajectoryData
		{
			std::string subsystem;
			std::vector<std::complex<double>> StateTrace;

			TrajectoryData(std::string name)
				:subsystem(name)
			{
				StateTrace = {};
			}
		};

		std::vector<std::pair<double, std::vector<TrajectoryData>>> trajectory;
		trajectory.push_back({0.0, {}});
		for (auto i = SubSystemSpins.cbegin(); i != SubSystemSpins.cend(); i++)
		{
			trajectory[0].second.push_back(TrajectoryData(i->first));
			// Get the superspace result vector and convert it back to the native Hilbert space
			arma::cx_mat rho_result;
			arma::cx_vec rho_result_vec;
			rho_result_vec = rho0.rows(nextDimension, nextDimension + SpinSpace->second->SpaceDimensions() - 1);
			//std::cout << rho_result_vec << std::endl;
			if (!SpinSpace->second->OperatorFromSuperspace(rho_result_vec, rho_result))
			{
				//this->Log() << "ERROR: Failed to convert resulting superspace-vector back to native Hilbert space for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}

			// Get the results
			this->GatherResults(rho_result, *(SpinSpace->first), *(SpinSpace->second), trajectory[0].second[i - SubSystemSpins.begin()].StateTrace);

			// Move on to next spin space
			nextDimension += SpinSpace->second->SpaceDimensions();
		}
		this->Data() << std::endl;

		// We need the propagator
		this->Log() << "Calculating the propagator..." << std::endl;
		//this->timestep = 1e-4;
		arma::cx_mat P = arma::expmat(arma::conv_to<arma::cx_mat>::from(L) * this->timestep);
		//arma::cx_mat P = arma::conv_to<arma::cx_mat>::from(L) * this->timestep;
	
		// Perform the calculation
		this->Log() << "Ready to perform calculation." << std::endl;
		unsigned int steps = static_cast<unsigned int>(std::abs(this->totaltime / this->timestep));
		for (unsigned int n = 1; n <= steps; n++)
		{
			// Write first part of the data output
			this->Data() << this->RunSettings()->CurrentStep() << " ";
			double time = static_cast<double>(n) * this->timestep;
			this->Data() << time << " ";
			this->WriteStandardOutput(this->Data());
			trajectory.push_back({time, {}});

			// Propagate (use special scope to be able to dispose of the temporary vector asap)
			{
				arma::cx_vec tmp = P * rho0;
				rho0 = tmp;
			}
			// Retrieve the resulting density matrix for each spin system and output the results
			nextDimension = 0;
			for (auto i = SubSystemSpins.cbegin(); i != SubSystemSpins.cend(); i++)
			{
				// Get the superspace result vector and convert it back to the native Hilbert space
				trajectory[n].second.push_back(TrajectoryData(i->first));
				arma::cx_mat rho_result;
				arma::cx_vec rho_result_vec;
				rho_result_vec = rho0.rows(nextDimension, nextDimension + SpinSpace->second->SpaceDimensions() - 1);
				if (!SpinSpace->second->OperatorFromSuperspace(rho_result_vec, rho_result))
				{
					//this->Log() << "ERROR: Failed to convert resulting superspace-vector back to native Hilbert space for spin system \"" << i->first->Name() << "\"!" << std::endl;
					return false;
				}

				// Get the results
				this->GatherResults(rho_result, *(SpinSpace->first), *(SpinSpace->second), trajectory[n].second[i - SubSystemSpins.begin()].StateTrace);

				// Move on to next spin space
				nextDimension += SpinSpace->second->SpaceDimensions();
			}

			// Terminate the line in the data file after iteration through all spin systems
			this->Data() << std::endl;
		}

		this->Log() << "Done with calculation." << std::endl;
		this->Log() << "Calculating Total Yield." << std::endl;
		//this->Data() << "\n" << "Yeilds ";
		std::vector<double> time;
		std::vector<double> yields;

		for(int i = 0; i < SpinSpace->first->States().size(); i++)
		{
			yields.push_back(0.0);
		}

		for(int i = 0; i < SubSystemSpins.size(); i++)
		{
			std::vector<std::vector<std::complex<double>>> data;
			for(auto e = trajectory.begin(); e != trajectory.end(); e++)
			{
				if(i == 0)
				{
					time.push_back(e->first);
				}

				for(auto a = e->second.begin(); a != e->second.end(); a++)
				{
					if(a->subsystem != SubSystemSpins[i].first)
					{
						continue;
					}

					data.push_back(a->StateTrace);
					break;
				}
			}

			for(auto e = SubSystemsTransitions.begin(); e != SubSystemsTransitions.end(); e++)
			{
				if(e->type != 0)
				{
					continue;
				}
				if(e->source != SubSystemSpins[i].first)
				{
					continue;
				}
				int index = 0;
				for(auto s : SpinSpace->first->States())
				{
					if(e->transition->SourceState() != s)
					{
						index++;
						continue;
					}
					break;
				}
				double rate = e->transition->Rate();
				double yield = 0;
				std::vector<std::complex<double>> state_data;
				for(auto a : data)
				{
					state_data.push_back(a[index]);
				}
				StateYield(rate, yield, state_data, time);
				yields[index] += yield;
			}
		}
		for(auto y : yields)
		{
			this->Data() << y << " ";
		}
		this->Data() << "\n" << std::endl;

		return true;
	}

	// Gathers and outputs the results from a given time-integrated density operator
	void TaskMultiRadicalPairSSTimeEvo::GatherResults(const arma::cx_mat &_rho, const SpinAPI::SpinSystem &_system, const SpinAPI::SpinSpace &_space, std::vector<std::complex<double>>& traj)
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
			//std::cout << P << std::endl;
			double tr = std::abs(arma::trace(P * _rho));
			this->Data() << tr << " ";
			//this->Data() << tr.real() << " ";
			traj.push_back(tr);
		}
	}

	void TaskMultiRadicalPairSSTimeEvo::StateYield(double _rate, double& _yeild, const std::vector<std::complex<double>>& _traj, std::vector<double>& _time)
	{
		auto f = [](double frac, double t, double kr) {return frac * std::exp(-kr * t); };
		std::vector<double> ylist;
		for(int i = 0; i < _traj.size(); i++)
		{
			//ylist.push_back(f(_traj[i].real(), _time[i], _rate));
			ylist.push_back(_rate * _traj[i].real());
		}
		_yeild = _rate * simpson_integration(_time, ylist);
	}

	double TaskMultiRadicalPairSSTimeEvo::simpson_integration(std::vector<double> x_list, std::vector<double> y_list)
	{
		double area = 0;
		for (int i = 0; i < x_list.size()-1; i++)
		{
			double diff = x_list[i + 1] - x_list[i];
			double ab = y_list[i] + y_list[i + 1];

			area = area + (ab * 0.5) * diff;
		}
		return area;
	}

	bool TaskMultiRadicalPairSSTimeEvo::GenerateHamiltonian(const std::vector<SpinAPI::interaction_ptr> interactions, arma::sp_cx_mat& H, int dimension, std::shared_ptr<SpinAPI::SpinSpace> SpinSystem)
	{
		if(interactions.size() < 1)
		{
			H = arma::sp_cx_mat(dimension, dimension);
			return true;
		}

		// Get the first interaction contribution
		auto i = interactions.cbegin();
		arma::sp_cx_mat tmp;
		arma::sp_cx_mat result; 
		if (!SpinSystem->InteractionOperator((*i), result))
			return false;

		// We have already used the first interaction
		i++;
		// Loop through the rest
		for (; i != interactions.cend(); i++)
		{
			// Attempt to get the matrix representing the Interaction object in the spin space
			if (!SpinSystem->InteractionOperator((*i), tmp))
				return false;
			result += tmp;
			//tmp.print();
		}

		H = result;
		//std::cout << H << std::endl;
		return true;
	}

	bool TaskMultiRadicalPairSSTimeEvo::GenerateReactionOperator(const std::vector<SubSystemTransition> transitions, arma::sp_cx_mat& K, int dimension, std::shared_ptr<SpinAPI::SpinSpace> SpinSystem, std::string source)
	{
		if(transitions.size() < 1)
		{
			K = arma::sp_cx_mat(dimension,dimension);
			return true;
		}

		// Get the first transition contribution
		auto i = transitions.cbegin();
		arma::sp_cx_mat tmp;
		arma::sp_cx_mat result;
		if(i->source == source)
		{
			if (!SpinSystem->ReactionOperator(i->transition, result, SpinAPI::ReactionOperatorType::Unspecified))
				return false;
		}

		//result.print();

		// We have already used the first transition
		i++;

		// Loop through the rest
		for (; i != transitions.cend(); i++)
		{
			if(i->source != source)
			{
				continue;
			}
			// Attempt to get the matrix representing the reaction operator in the spin space
			if (!SpinSystem->ReactionOperator(i->transition, tmp, SpinAPI::ReactionOperatorType::Unspecified))
				return false;
			result += tmp;
			//tmp.print();
		}

		K = result;
		//K.print();
		return true;

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

		ValidateSubSystems();

		return true;
	}

	bool TaskMultiRadicalPairSSTimeEvo::ValidateSubSystems()
	{
		int SubSystemsNum = 0;
		if(!this->Properties()->Get("subsystems", SubSystemsNum))
		{
			this->Log() << "Error: Number of sub systems not specified. Correct syntax is e.g. subsystems: 2" << std::endl;
			return false;
		}

		std::string SubSystemNames; 
		if(!this->Properties()->Get("subsystemnames", SubSystemNames))
		{
			this->Log() << "Note: subsystemnames not specified. Asssuming names are system1, system2, system3...." << std::endl;
		}

		std::vector<std::string> SubSystems;
		{
			std::stringstream ss(SubSystemNames);
			while(ss.good())
			{
				std::string Name;
				std::getline(ss,Name, ',');
				std::transform(Name.begin(), Name.end(), Name.begin(), [](unsigned char c){return std::tolower(c);});
				SubSystems.push_back(Name);
			}
		}

		//Verify that the subsytems are valid by checking if all the objects have been loaded
		int num = 0;
		for(auto system : SubSystems)
		{
			std::vector<std::string> SubSystemTemp;
			std::string SubSytemDefinition = "";
			if(!this->Properties().get()->Get(system, SubSytemDefinition))
			{
				this->Log() << "Error: Failed to Load SubSystem " << system << std::endl;
				std::cout << "Error: Failed to Load SubSystem " << system << std::endl;
				num++;
				continue;
			}
		}
		return true;
	}
    
	bool TaskMultiRadicalPairSSTimeEvo::RungeKutta4(arma::sp_cx_mat &L, arma::cx_vec &RhoNaught, arma::cx_vec &drhodt, double timestep)
    {
        return false;
    }
    
	arma::cx_vec TaskMultiRadicalPairSSTimeEvo::ComputeRhoDot(arma::sp_cx_mat &L, arma::cx_vec &K, amra::cx_vec &RhoNaugt)
    {
        return arma::cx_vec();
    }
}