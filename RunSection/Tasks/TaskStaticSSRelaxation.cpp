/////////////////////////////////////////////////////////////////////////
// TaskStaticSSRelaxation implementation (RunSection module)
// 
// (c) 2018 by Claus N.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TaskStaticSSRelaxation.h"
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "Spin.h"
#include "ObjectParser.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskStaticSSRelaxation Constructors and Destructor
	// -----------------------------------------------------
	TaskStaticSSRelaxation::TaskStaticSSRelaxation(const MSDParser::ObjectParser& _parser, const RunSection& _runsection)	: BasicTask(_parser,_runsection), relaxationRate(0.0), relaxingSpins(), productYieldsOnly(false)
	{
		
	}
	
	TaskStaticSSRelaxation::~TaskStaticSSRelaxation()
	{
		
	}
	// -----------------------------------------------------
	// TaskStaticSSRelaxation protected methods
	// -----------------------------------------------------
	bool TaskStaticSSRelaxation::RunLocal()
	{
		this->Log() << "Running method StaticSS." << std::endl;
		
		// If this is the first step, write first part of header to the data file
		if(this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}
		
		// Loop through all SpinSystems
		auto systems = this->SpinSystems();
		for(auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Make sure we have an initial state
			auto initial_states = (*i)->InitialState();
			if(initial_states.size() < 1)
			{
				this->Log() << "Skipping SpinSystem \"" << (*i)->Name() << "\" as no initial state was specified." << std::endl;
				continue;
			}
			
			this->Log() << "\nStarting with SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
			
			// Obtain a SpinSpace to describe the system
			SpinAPI::SpinSpace space(*(*i));
			arma::cx_mat rho0;
			space.UseSuperoperatorSpace(true);
			
			// Get the initial state
			for(auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
			{
				arma::cx_mat tmp_rho0;
				if(!space.GetState(*j, tmp_rho0))
				{
					this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
					continue;
				}
				if(j == initial_states.cbegin())
					rho0 = tmp_rho0;
				else
					rho0 += tmp_rho0;
			}
			rho0 /= arma::trace(rho0);	// The density operator should have a trace of 1
			
			// Convert initial state to superoperator space
			arma::cx_vec rho0vec;
			if(!space.OperatorToSuperspace(rho0, rho0vec))
			{
				this->Log() << "Failed to convert initial state density operator to superspace." << std::endl;
				continue;
			}
			
			// Get the Hamiltonian
			arma::sp_cx_mat H;
			if(!space.Hamiltonian(H))
			{
				this->Log() << "Failed to obtain Hamiltonian in superspace." << std::endl;
				continue;
			}
			
			// Get a matrix to collect all the terms (the total Liouvillian)
			arma::sp_cx_mat A = arma::cx_double(0.0, -1.0) * H;
			
			// Get the reaction operators, and add them to "A"
			arma::sp_cx_mat K;
			if(!space.TotalReactionOperator(K))
			{
				this->Log() << "Warning: Failed to obtain matrix representation of the reaction operators!" << std::endl;
			}
			A -= K;
			
			// ----------------------------------------------------------------
			// Setup relaxation operator
			for(auto j = this->relaxingSpins.cbegin(); j != this->relaxingSpins.cend(); j++)
			{
				// Attempt to find the spin in the spin system
				auto spin = (*i)->spins_find( *j );
				if(spin == nullptr)
					continue;
				this->Log() << "Found spin \"" << (*j) << "\" in spin system \"" << (*i)->Name() << "\"." << std::endl;
				
				// Obtain spin operators
				arma::sp_cx_mat Sx;
				arma::sp_cx_mat Sy;
				arma::sp_cx_mat Sz;
				if(!space.CreateOperator(spin->Sx(), spin, Sx) || !space.CreateOperator(spin->Sy(), spin, Sy) || !space.CreateOperator(spin->Sz(), spin, Sz))
				{
					this->Log() << "Failed to obtain spin operators for spin \"" << (*j) << "\"." << std::endl;
					continue;
				}
				
				// Convert spin operators to superoperator space
				arma::sp_cx_mat SSx1;
				arma::sp_cx_mat SSy1;
				arma::sp_cx_mat SSz1;
				arma::sp_cx_mat SSx2;
				arma::sp_cx_mat SSy2;
				arma::sp_cx_mat SSz2;
				arma::sp_cx_mat SSx3;
				arma::sp_cx_mat SSy3;
				arma::sp_cx_mat SSz3;
				if(!space.SuperoperatorFromOperators(Sx, Sx, SSx1) || !space.SuperoperatorFromOperators(Sy, Sy, SSy1) || !space.SuperoperatorFromOperators(Sz, Sz, SSz1)
					|| !space.SuperoperatorFromLeftOperator(Sx*Sx, SSx2) || !space.SuperoperatorFromLeftOperator(Sy*Sy, SSy2) || !space.SuperoperatorFromLeftOperator(Sz*Sz, SSz2)
					|| !space.SuperoperatorFromRightOperator(Sx*Sx, SSx3) || !space.SuperoperatorFromRightOperator(Sy*Sy, SSy3) || !space.SuperoperatorFromRightOperator(Sz*Sz, SSz3))
				{
					this->Log() << "Failed to obtain spin operators for spin \"" << (*j) << "\"." << std::endl;
					continue;
				}
				
				// Add the relaxation contribution to the Liouvillian
				A += 0.5 * this->relaxationRate * (2*(SSx1 + SSy1 + SSz1) - (SSx2 + SSy2 + SSz2) - (SSx3 + SSy3 + SSz3));
			}
			// ----------------------------------------------------------------
			
			// Perform the calculation
			this->Log() << "Ready to perform calculation." << std::endl;
			arma::cx_vec result = solve(arma::conv_to<arma::cx_mat>::from(A), rho0vec);
			this->Log() << "Done with calculation." << std::endl;
			
			// Convert the resulting density operator back to its Hilbert space representation
			if(!space.OperatorFromSuperspace(result, rho0))
			{
				this->Log() << "Failed to convert resulting superspace-vector back to native Hilbert space." << std::endl;
				continue;
			}
			
			// Obtain the results
			arma::cx_mat P;
			this->Data() << this->RunSettings()->CurrentStep() << " ";
			this->WriteStandardOutput(this->Data());
			
			// There are two result modes - either write results per transition or for each defined state
			if(this->productYieldsOnly)
			{
				// Loop through all defind transitions
				auto transitions = (*i)->Transitions();
				for(auto j = transitions.cbegin(); j != transitions.cend(); j++)
				{
					// Make sure that there is a state object
					if((*j)->SourceState() == nullptr)
						continue;
					
					if(!space.GetState((*j)->SourceState(), P))
					{
						this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
						continue;
					}
					
					// Return the yield for this transition
					this->Data() << (*j)->Rate() * std::abs(arma::trace(P * rho0)) << " ";
				}
			}
			else
			{
				// Loop through all states
				auto states = (*i)->States();
				for(auto j = states.cbegin(); j != states.cend(); j++)
				{
					if(!space.GetState((*j), P))
					{
						this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
						continue;
					}
					
					// Return the yield for this state - note that no reaction rates are included here.
					this->Data() << std::abs(arma::trace(P * rho0)) << " ";
				}
			}
			
			this->Log() << "\nDone with SpinSystem \"" << (*i)->Name() << "\"" << std::endl;
		}
		
		// Terminate the line in the data file after iteration through all spin systems
		this->Data() << std::endl;
		
		return true;
	}
	
	// Writes the header of the data file (but can also be passed to other streams)
	void TaskStaticSSRelaxation::WriteHeader(std::ostream& _stream)
	{
		_stream << "Step ";
		this->WriteStandardOutputHeader(_stream);
		
		// Get header for each spin system
		auto systems = this->SpinSystems();
		for(auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Should yields be written per transition or per defined state?
			if(this->productYieldsOnly)
			{
				// Write each transition name
				auto transitions = (*i)->Transitions();
				for(auto j = transitions.cbegin(); j != transitions.cend(); j++)
					_stream << (*i)->Name() << "." << (*j)->Name() << ".yield ";
			}
			else
			{
				// Write each state name
				auto states = (*i)->States();
				for(auto j = states.cbegin(); j != states.cend(); j++)
					_stream << (*i)->Name() << "." << (*j)->Name() << " ";
			}
		}
		_stream << std::endl;
	}
	
	// Validation
	bool TaskStaticSSRelaxation::Validate()
	{
		this->Properties()->Get("transitionyields",this->productYieldsOnly);
		
		this->Properties()->Get("relaxationrate",this->relaxationRate);
		
		this->Properties()->GetList("relaxingspins",this->relaxingSpins);
		
		return true;
	}
	// -----------------------------------------------------
}

