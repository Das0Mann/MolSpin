/////////////////////////////////////////////////////////////////////////
// State class (SpinAPI Module)
// ------------------
// The State class represents a spin state.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_State
#define MOD_SpinAPI_State

#include <memory>
#include <armadillo>
#include <unordered_map>
#include "SpinAPIfwd.h"
#include "SpinSystem.h"
#include "Function.h"

namespace SpinAPI
{
	class State
	{
		// Aliases to improve readability
		using StateSeries = std::vector<std::pair<int, arma::cx_double>>; // A collection of "mz" values and state factors.

		using StatePair = std::pair<spin_ptr, StateSeries>; // A spin_ptr and the corresponding StateSeries.
															// For non-entangled states the vector should contain 1 element.

		using CompleteState = std::vector<StatePair>; // A collection of states that are entangled together,
													  // or only 1 state if it is not entangled with any other states.
	private:
		// Data members
		std::shared_ptr<MSDParser::ObjectParser> properties; // Use a pointer to the object to minimize compilation dependencies
		std::vector<CompleteState> substates;				 // A list of complete states, i.e. see the alias definition above
		std::vector<std::shared_ptr<Function>> Functions; 	 // A list of functions that act as factors before states, e.g. cos or sin
		std::vector<int> BracketDepth;						 // A description of the state when it comes to how to apply factors
		std::vector<arma::cx_double> InitialFactors; 				 // stores the initial factors (see UpdateFactors for reasoning)
		std::unordered_map<std::string, double> Variables; 	// Variable name and function and it's value(stored as a void* for type casting)
		bool isValid;

		// Private methods
		bool ParseMz(const std::string &, int &);							  // Gets the mz value from the string
		bool ParseFactor(const std::string &, arma::cx_double &);			  // Converts string into complex number (cx_double) if possible, i.e. the factor for the SpinState
		bool ParseStates(const std::vector<spin_ptr> &, const std::string &); // Parses the state string to fill the substates vector

		const CompleteState *FindState(const spin_ptr &, StateSeries **_state = nullptr) const; // Returns a CompleteState pointer if the spin_ptr was found, nullptr otherwise
																								// Sets the StateSeries pointer if the pointer-to-pointer is not nullptr

		// Private methods to create ActionTargets
		std::vector<RunSection::NamedActionScalar> CreateActionScalars(const std::string &); 

		bool UpdateFactors();   

	public:
		// Constructors / Destructors
		State(std::string, std::string); // Normal constructor
		State(const State &);			 // Copy-constructor
		~State();						 // Destructor

		// Operators
		const State &operator=(const State &); // Copy-assignment

		// Name and validation
		std::string Name() const;
		bool IsValid() const;

		// Allow access to custom properties to be used for custom tasks
		std::shared_ptr<const MSDParser::ObjectParser> Properties() const;

		// Attempts to parse the State object using spins from the given SpinSystem, and sets this->system if successful
		// Cannot be called it a system was already assigned successfully (returns false in that case)
		bool ParseFromSystem(const SpinSystem &);

		// Public methods
		bool IsCoupled(const spin_ptr &) const;								// Checks whether the spin is in a coupled/entangled state
		bool IsComplete(const std::vector<spin_ptr> &) const;				// Checks whether any spins in the vector are coupled/entangled to any spins outside the vector
		bool CompleteSet(std::vector<spin_ptr> &) const;					// Extends the set such that all spins inside the set can only be coupled/entangled to other spins inside the set
		bool GetStateSeries(const spin_ptr &, StateSeries &) const;			// Sets the StateSeries for the given spin if it is well-defined (not entangled) for the state, returns false otherwise
		bool GetStateSeries(const spin_ptr &, StateSeries &, bool &) const; // Sets the bool to true if the spin was coupled, otherwise same as above
		bool GetCompleteState(const spin_ptr &, CompleteState &) const;		// Sets the CompleteState for the given spin if it is defined for the state, returns false otherwise

		// The print function writes the contents of the State object to an output-stream
		// This function is intended for testing/debugging purposes and is specific to the class implementation
		void Print(std::ostream &, unsigned int _name_width = 20, unsigned int _data_width = 8) const;

		// Public method for creating ActionTargets
		void GetActionTargets(std::vector<RunSection::NamedActionScalar> &, std::vector<RunSection::NamedActionVector> &, const std::string &);

		//Applies all the functions to the factors
		bool Update();
	
	};

	// Define alias for state-pointers
	using state_ptr = std::shared_ptr<State>;

	bool CheckActionScalarVariable(const double &);


}

#endif
