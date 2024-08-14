/////////////////////////////////////////////////////////////////////////
// State implementation (SpinAPI Module)
// ------------------
// The State class represents a spin state.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iomanip>
#include "ObjectParser.h"
#include "Spin.h"
#include "State.h"
#include "SpinSystem.h"

namespace SpinAPI
{
	// -----------------------------------------------------
	// State Constructors and Destructor
	// -----------------------------------------------------
	State::State(std::string _name, std::string _contents) : properties(std::make_shared<MSDParser::ObjectParser>(_name, _contents)), substates(), isValid(false)
	{
	}

	State::State(const State &_state) : properties(std::make_shared<MSDParser::ObjectParser>(*(this->properties))), substates(), isValid(_state.isValid)
	{
	}

	State::~State()
	{
	}
	// -----------------------------------------------------
	// Operators
	// -----------------------------------------------------
	// Copy-assignment
	const State &State::operator=(const State &_state)
	{
		this->properties = std::make_shared<MSDParser::ObjectParser>(*(this->properties));
		this->substates = _state.substates;

		return (*this);
	}
	// -----------------------------------------------------
	// Methods to convert strings to numbers
	// -----------------------------------------------------
	// Parse integer from string
	// Note that "mz" is stored in units of "1/2"
	bool State::ParseMz(const std::string &_mz, int &_out)
	{
		try
		{
			// Check whether mz is written in the form "1/2", "2/2", "3/2", etc.
			if (_mz.size() > 2 && (*(_mz.cend() - 1) == '2' && *(_mz.cend() - 2) == '/'))
			{
				_out = std::stoi(_mz.substr(0, _mz.size() - 2).c_str());
			}
			else
			{
				// Use a factor of two if mz is not specified with "/2"
				_out = 2 * std::stoi(_mz.c_str());
			}
		}
		catch (const std::exception &)
		{
			return false;
		} // In case of exceptions from std::stoi

		return true;
	}

	// Parse complex number from string
	bool State::ParseFactor(const std::string &_factor, arma::cx_double &_out)
	{
		// An empty string (no explicit factor) is just a factor of 1
		// Also check for other similar simple cases
		if (_factor.empty() || _factor.compare("+") == 0)
		{
			_out = arma::cx_double(1, 0);
			return true;
		}
		else if (_factor.compare("-") == 0)
		{
			_out = arma::cx_double(-1, 0);
			return true;
		}
		else if (_factor.compare("+i") == 0 || _factor.compare("i") == 0)
		{
			_out = arma::cx_double(0, 1);
			return true;
		}
		else if (_factor.compare("-i") == 0)
		{
			_out = arma::cx_double(0, -1);
			return true;
		}

		double re = 0.0;
		double im = 0.0;
		std::string buffer = "";
		bool isImaginary = false;
		bool foundReal = false;
		bool foundImag = false;

		try
		{
			// Iterate through the string
			for (auto i = _factor.cbegin(); i != _factor.cend(); i++)
			{
				// Ignore parantheses here
				if ((*i) == '(' || (*i) == ')')
					continue;

				if ((*i) == '+' || (*i) == '-')
				{
					// Do we have a number, or is this the start of a number?
					if (!buffer.empty())
					{
						// Save the number we have
						if (isImaginary)
						{
							im = std::stod(buffer);
							foundImag = true; // This bool is tested when checking for character 'i', so we don't need to consider if an imaginary part was found previously
						}
						else
						{
							// There should be only a single contribution to the real part
							if (foundReal)
								return false;

							re = std::stod(buffer);
							foundReal = true;
						}

						// Reset
						buffer = "";
						isImaginary = false;
					}

					buffer += (*i);
				}
				else if ((*i) == 'i' || (*i) == 'I')
				{
					// There should be only one 'i' per number, and only a single imaginary part specified
					if (isImaginary || foundImag)
						return false;

					isImaginary = true; // TODO: Check that 'i' isn't in the middle of the number
				}
				else
				{
					buffer += (*i);
				}
			}

			// Do we have a number
			if (!buffer.empty())
			{
				// Save the number we have
				if (isImaginary)
				{
					im = std::stod(buffer);
					foundImag = true; // This bool is tested when checking for character 'i', so we don't need to consider if an imaginary part was found previously
				}
				else
				{
					// There should be only a single contribution to the real part
					if (foundReal)
						return false;

					re = std::stod(buffer);
					foundReal = true;
				}
			}
		}
		catch (const std::exception &)
		{
			return false;
		} // Catch exceptions from the uses of std::stod

		// Check whether anything was found
		if (!foundReal && !foundImag)
			return false;

		// Return the result
		_out = arma::cx_double(re, im);
		return true;
	}
	// -----------------------------------------------------
	// Private parse function to get SpinState structs
	// -----------------------------------------------------
	// Parse states from a string for the given spins
	bool State::ParseStates(const std::vector<spin_ptr> &_spins, const std::string &_states)
	{
		CompleteState newState;

		// Prepare newState and make sure that none of the spin are already accounted for
		{
			StateSeries emptySeries;

			for (auto i = _spins.cbegin(); i != _spins.cend(); i++)
			{
				// Check that the same spin is not specified twice within the same state
				if (std::find(i + 1, _spins.cend(), (*i)) != _spins.cend())
					return false;

				// Check whether the state was already specified for the spin
				auto tmpState = this->FindState((*i));
				if (tmpState != nullptr)
					return false;

				// Create StatePair by copying empty StateSeries object
				newState.push_back(StatePair((*i), emptySeries));
			}
		}

		std::string buffer = "";
		arma::cx_double factor;
		arma::cx_double PreFactor = 0;
		int mz = 0;
		bool inState = false;
		auto currentSpinPair = newState.begin();
		bool brackets = false;

		// Loop through all characters in the string
		for (auto i = _states.cbegin(); i != _states.cend(); i++)
		{
			if (!inState && (*i) == '|')
			{
				// Attempt to parse the factor
				if (!this->ParseFactor(buffer, factor))
					return false;

				if(brackets)
				{
					factor = factor * PreFactor;
				}

				// We are now inside a state (not the factor), so reset the buffer
				inState = true;
				buffer = "";
			}
			else if (inState && (*i) == ',')
			{
				// Make sure that we still have a spin left to assign an "mz" value to
				if (currentSpinPair == newState.end())
					return false;

				// Attempt to parse the "mz" value
				if (!ParseMz(buffer, mz))
					return false;

				// Check whether the mz value is allowed
				if (mz > currentSpinPair->first->S() || mz < -currentSpinPair->first->S())
				{
					std::cout << "ERROR: Value " << mz << "/2 for mz is not allowed for a spin of " << currentSpinPair->first->S() << "/2!" << std::endl;
					return false;
				}

				// Extend the StateSeries with a new pair of "mz" and "factor" values
				currentSpinPair->second.push_back(std::pair<int, arma::cx_double>(mz, factor));

				// Reset buffer and prepare reading next mz value
				buffer = "";
				++currentSpinPair;
			}
			else if (inState && (*i) == '>')
			{
				// Make sure that we still have a spin left to assign an "mz" value to
				if (currentSpinPair == newState.end())
					return false;

				// Attempt to parse the "mz" value
				if (!ParseMz(buffer, mz))
					return false;

				// Check whether the mz value is allowed
				if (mz > currentSpinPair->first->S() || mz < -currentSpinPair->first->S())
				{
					std::cout << "ERROR: Value " << mz << "/2 for mz is not allowed for a spin of " << currentSpinPair->first->S() << "/2!" << std::endl;
					return false;
				}

				// Extend the StateSeries with a new pair of "mz" and "factor" values
				currentSpinPair->second.push_back(std::pair<int, arma::cx_double>(mz, factor));

				// Reset buffer and prepare to read the next state
				buffer = "";
				inState = false;
				++currentSpinPair;

				// Make sure that we used exactly all the spins
				if (currentSpinPair != newState.end())
					return false;

				// Reset the CompleteState iterator
				currentSpinPair = newState.begin();
			}
			else if (!inState && (*i) == '(')
			{
				if (!this->ParseFactor(buffer, factor))
					return false;
				
				buffer = "";

				PreFactor = factor;
				brackets = true;
			}
			else if (!inState && (*i) == ')')
			{
				PreFactor = 0;
				brackets = false;
			}
			else
			{
				buffer += (*i);
			}
		}

		this->substates.push_back(newState);
		return true;
	}
	// -----------------------------------------------------
	// Public method that gets information from the properties
	// object and prepares it for a call to ParseStates
	// -----------------------------------------------------
	bool State::ParseFromSystem(const SpinSystem &_system)
	{
		// Clean-up if we attempted to create the state previously
		if (!this->substates.empty())
			this->substates.clear();

		// Get a list of coupled spin state specifications
		auto coupled = this->properties->GetFunction("spins");

		// Loop through the list
		for (auto i = coupled.cbegin(); i != coupled.cend(); i++)
		{
			// Helper variables
			std::vector<spin_ptr> spinlist;
			std::string buffer = "";
			spin_ptr tempSpin = nullptr;

			// Loop through the string
			for (auto j = i->first.cbegin(); j != i->first.cend(); j++)
			{
				// The string is comma-separated
				if ((*j) == ',')
				{
					// Get the spin with the requested name, or abort if it fails
					tempSpin = _system.spins_find(buffer);
					if (tempSpin == nullptr)
						return false;

					// Save it to the list of found spins
					spinlist.push_back(tempSpin);

					// Reset the buffer
					buffer = "";
				}
				else
				{
					// If not a comma, it must be part of the spin name
					buffer += (*j);
				}
			}

			// Don't forget to also get the last spin if the buffer is not empty
			if (!buffer.empty())
			{
				tempSpin = _system.spins_find(buffer);
				if (tempSpin == nullptr)
					return false;
				spinlist.push_back(tempSpin);
			}

			// We are done with the string, so call the method that parses the state itself
			// Abort if the state could not be parsed
			if (!this->ParseStates(spinlist, i->second))
				return false;
		}

		// Single-spin states are simpler to parse
		auto uncoupled = this->properties->GetFunction("spin");
		for (auto i = uncoupled.cbegin(); i != uncoupled.cend(); i++)
		{
			// Prepare a list
			std::vector<spin_ptr> spinlist;

			// Get the spin with the requested name, or abort if it fails
			spin_ptr tempSpin = _system.spins_find(i->first);
			if (tempSpin == nullptr)
				return false;

			// Save it to the list of found spins, and call the parse function
			spinlist.push_back(tempSpin);

			// Abort if the state could not be parsed
			if (!this->ParseStates(spinlist, i->second))
				return false;
		}

		// Don't waste any unnecessary memory
		this->substates.shrink_to_fit();

		// If we were successful, the state is now valid
		this->isValid = true;
		return true;
	}
	// -----------------------------------------------------
	// Private methods to search for internal state objects
	// -----------------------------------------------------
	// Returns a CompleteState pointer if the spin_ptr was found, nullptr otherwise
	// Sets the StateSeries pointer if the pointer-to-pointer is not nullptr
	const State::CompleteState *State::FindState(const spin_ptr &_spin, StateSeries **_state) const
	{
		// Loop through each CompleteState
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			// Loop through each StatePair within the CompleteState
			for (auto j = i->cbegin(); j != i->cend(); j++)
			{
				// Check if we found the right StatePair
				if (j->first == _spin)
				{
					// If we have a pointer to a StateSeries pointer, set it
					if (_state != nullptr)
						(*_state) = const_cast<StateSeries *>(&(j->second)); // TODO: Get rid of const_cast

					// Return a reference to the CompleteState
					return &(*i);
				}
			}
		}

		return nullptr;
	}
	// -----------------------------------------------------
	// Name and validation
	// -----------------------------------------------------
	std::string State::Name() const
	{
		return this->properties->Name();
	}

	bool State::IsValid() const
	{
		return this->isValid;
	}
	// -----------------------------------------------------
	// Access to custom properties
	// -----------------------------------------------------
	std::shared_ptr<const MSDParser::ObjectParser> State::Properties() const
	{
		return this->properties;
	}
	// -----------------------------------------------------
	// Public Get State methods
	// -----------------------------------------------------
	// Checks whether the spin is coupled to other spins (entanglement)
	bool State::IsCoupled(const spin_ptr &_spin) const
	{
		// Loop through each CompleteState
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			// There has to be at least two spins in the CompleteState for them to be coupled/entangled
			if ((*i).size() < 2)
			{
				continue;
			}

			// Loop through each StatePair within the CompleteState
			for (auto j = i->cbegin(); j != i->cend(); j++)
			{
				// Check if we found the right spin
				if (j->first == _spin)
				{
					// The spin is coupled/entangled if the CompleteState containing the spin also contains other spins
					if ((*i).size() > 1)
						return true;
				}
			}
		}

		return false;
	}

	// Checks whether any spins in the vector are coupled/entangled to any spins outside the vector
	bool State::IsComplete(const std::vector<spin_ptr> &_spins) const
	{
		unsigned int spinsFound = 0;

		// Loop through each CompleteState
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			// There has to be at least two spins in the CompleteState for them to be coupled/entangled
			if ((*i).size() < 2)
			{
				continue;
			}

			spinsFound = 0; // Reset counter

			// Loop through each StatePair within the CompleteState
			for (auto j = i->cbegin(); j != i->cend(); j++)
			{
				// Is the current spin within the collection?
				if (std::find(_spins.cbegin(), _spins.cend(), j->first) != _spins.cend())
					++spinsFound;
			}

			// Did the CompleteState contain some but not all spins in the _spins collection?
			// If so, the _spins collection is not complete!
			if (spinsFound > 0 && spinsFound < (*i).size())
				return false;
		}

		return true;
	}

	// Extends the set such that all spins inside the set can only be coupled/entangled to other spins inside the set
	bool State::CompleteSet(std::vector<spin_ptr> &_spins) const
	{
		// List of entangled spins
		std::vector<spin_ptr> entangledspins;
		entangledspins.reserve(_spins.size()); // Only allocate the same size as the spin list - this could be enough since only entangled spins will be added

		// Loop through each CompleteState
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			// There has to be at least two spins in the CompleteState for them to be coupled/entangled
			if ((*i).size() < 2)
			{
				continue;
			}

			// Temporary/helper variables
			bool addCompleteStateSpins = false;
			std::vector<spin_ptr> templist;
			templist.reserve(i->size()); // We know the size, so allocate it from the beginning

			// Loop through each StatePair within the CompleteState
			for (auto j = i->cbegin(); j != i->cend(); j++)
			{
				// Get a list of all spins in the CompleteState
				templist.push_back(j->first);

				// Is the current spin within the collection?
				if (std::find(_spins.cbegin(), _spins.cend(), j->first) != _spins.cend())
					addCompleteStateSpins = true;
			}

			// If a spin in "_spins" is in the CompleteState, get all the spins in that CompleteState
			if (addCompleteStateSpins)
				entangledspins.insert(entangledspins.end(), templist.begin(), templist.end());
		}

		// Prepare for the insertion
		auto spincount = _spins.size(); // Get number of spins before changing the collection
		if (entangledspins.size() > _spins.capacity())
			_spins.reserve(entangledspins.size()); // Make sure that we don't have multiple allocations

		// Insert the missing spins; i.e. extend the spin subspace
		for (auto i = entangledspins.cbegin(); i != entangledspins.cend(); i++)
			if (std::find(_spins.cbegin(), _spins.cend(), (*i)) == _spins.cend())
				_spins.push_back((*i));

		// Return true if any spins were added
		return (_spins.size() > spincount);
	}

	// Gets the state series if the spin is NOT coupled/entangled
	// Returns false if the spin was coupled/entangled or if it was not found
	bool State::GetStateSeries(const spin_ptr &_spin, StateSeries &_state) const
	{
		bool isCoupled = false;
		return GetStateSeries(_spin, _state, isCoupled);
	}

	// Same as above, but also gives information about whether the spin is coupled/entangled
	bool State::GetStateSeries(const spin_ptr &_spin, StateSeries &_state, bool &_isCoupled) const
	{
		// Loop through each CompleteState
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			// Loop through each StatePair within the CompleteState
			for (auto j = i->cbegin(); j != i->cend(); j++)
			{
				// Check if we found the right spin
				if (j->first == _spin)
				{
					// Check whether the spin is coupled/entangled
					if ((*i).size() > 1)
					{
						_isCoupled = true;
						return false;
					}
					else
					{
						_isCoupled = false;
						_state = j->second;
						return true;
					}
				}
			}
		}

		_isCoupled = false;
		return false;
	}

	// Returns a CompleteState containing the spin
	bool State::GetCompleteState(const spin_ptr &_spin, CompleteState &_state) const
	{
		// Loop through each CompleteState
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			// Loop through each StatePair within the CompleteState
			for (auto j = i->cbegin(); j != i->cend(); j++)
			{
				// Check if we found the right spin
				if (j->first == _spin)
				{
					_state = (*i);
					return true;
				}
			}
		}

		return false;
	}

	// -----------------------------------------------------
	// Other public methods
	// -----------------------------------------------------
	// Prints the contents of the State object
	void State::Print(std::ostream &_stream, unsigned int _name_width, unsigned int _data_width) const
	{
		// Loop through CompleteState objects
		for (auto i = this->substates.cbegin(); i != this->substates.cend(); i++)
		{
			if (i->size() > 0)
			{
				// Write header with factors
				_stream << std::setw(_name_width) << this->properties->Name() << "|";
				for (auto j = (*i)[0].second.cbegin(); j != (*i)[0].second.cend(); j++)
					_stream << std::setw(_data_width) << j->second;
				_stream << std::endl;

				_stream << std::string(_name_width + _data_width * (*i)[0].second.size() + 1, '-') << std::endl;

				// Loop through StatePair objects
				for (auto j = i->cbegin(); j != i->cend(); j++)
				{
					_stream << std::setw(_name_width) << j->first->Name() << "|";

					// Loop through StateSeries object
					for (auto k = j->second.cbegin(); k != j->second.cend(); k++)
						_stream << std::setw(_data_width) << k->first;
					_stream << std::endl;
				}
			}
			else
			{
				_stream << "CompleteState object was empty!" << std::endl;
			}

			if (i + 1 != this->substates.cend())
				_stream << "\n";
		}
	}
	// -----------------------------------------------------
}
