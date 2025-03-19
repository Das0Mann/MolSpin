/////////////////////////////////////////////////////////////////////////
// SpinSpace class (SpinAPI Module)
// ------------------
// This source file contains methods for managing which spins and
// interactions that are handled by the SpinSpace class, i.e. the methods
// Add, Remove, Clear and Contains.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
namespace SpinAPI
{
	// -----------------------------------------------------
	// Spin Management: Add-methods
	// -----------------------------------------------------
	// Add spin if not already contained by the SpinSpace, and returns true if is was added
	bool SpinSpace::Add(const spin_ptr &_spin)
	{
		if (std::find(this->spins.cbegin(), this->spins.cend(), _spin) != this->spins.cend())
			return false;

		this->spins.push_back(_spin);
		return true;
	}

	// Adds all the spins from a given list, without allowing any duplicates. Returns true if any elements were added
	// NOTE: This method may change the order of the spins in the vector if and only if new elements were inserted.
	// NOTE furthermore that this may change the order of the basis vectors of the space! TODO: Fix this before it causes confusion.
	bool SpinSpace::Add(const std::vector<spin_ptr> &_list)
	{
		// Merge the two lists without keeping any duplicates
		std::vector<spin_ptr> tmpvec;										   // Make a temporary list...
		tmpvec.reserve(this->spins.size() + _list.size());					   // ...and reserve memory to hold the following three inserts
		tmpvec.insert(tmpvec.end(), _list.cbegin(), _list.cend());			   // Insert all the new elements from the list
		tmpvec.insert(tmpvec.end(), this->spins.cbegin(), this->spins.cend()); // Insert all the old elements
		std::sort(tmpvec.begin(), tmpvec.end());							   // Sort the elements
		auto newEnd = std::unique(tmpvec.begin(), tmpvec.end());			   // Remove duplicates (requires list to be sorted)
		tmpvec.resize(std::distance(tmpvec.begin(), newEnd));				   // Shrink vector size to get rid of removed objects (std::unique cannot do that)

		// Check if any spins were added (if the new collection contains more elements than the old)
		if (tmpvec.size() > this->spins.size())
		{
			this->spins = tmpvec;
			return true;
		}

		return false;
	}

	// Adds all the spins and interactions from the SpinSystem to the spin space.
	bool SpinSpace::Add(const std::shared_ptr<SpinSystem> &_system)
	{
		if (_system != nullptr)
			return (this->Add(_system->Spins()) & this->Add(_system->Interactions()) & this->Add(_system->Transitions()) & this->Add(_system->Pulses()));

		return false;
	}

	// Adds all the spins and interactions from the SpinSystem to the spin space.
	bool SpinSpace::Add(const SpinSystem &_system)
	{
		return (this->Add(_system.Spins()) & this->Add(_system.Interactions()) & this->Add(_system.Transitions()) & this->Add(_system.Pulses()));
	}

	// Add missing spins (if any) to make the space comeplete with respect to the state
	bool SpinSpace::CompleteSpace(const state_ptr &_state) // TODO: Implement
	{
		std::cout << "Method SpinSpace::CompleteSpace(const state_ptr&) is not implemented!" << std::endl;
		return false;
	}

	// Add missing spins (if any) to make the space comeplete with respect to the interaction
	bool SpinSpace::CompleteSpace(const interaction_ptr &_interaction)
	{
		return _interaction->CompleteSet(this->spins);
	}
	// -----------------------------------------------------
	// Spin Management: Remove-methods
	// -----------------------------------------------------
	// Remove single spin
	bool SpinSpace::Remove(const spin_ptr &_spin)
	{
		auto i = std::find(this->spins.cbegin(), this->spins.cend(), _spin);
		if (i != this->spins.cend())
		{
			this->spins.erase(i);
			return true;
		}

		return false;
	}

	// Remove collection of spins
	bool SpinSpace::Remove(const std::vector<spin_ptr> &_list) // TODO: Implement more efficient version of the method
	{
		bool removed_any = false;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			removed_any |= this->Remove(*i);
		return removed_any;
	}

	// Remove all spins from a given SpinSystem
	bool SpinSpace::Remove(const std::shared_ptr<SpinSystem> &_system)
	{
		if (_system == nullptr)
			return false;

		return this->Remove(_system->Spins());
	}

	// Same as above
	bool SpinSpace::Remove(const SpinSystem &_system)
	{
		return this->Remove(_system.Spins());
	}

	// Remove the given spin and all spins it is entangled to
	bool SpinSpace::RemoveSubspace(const spin_ptr &_spin, const state_ptr &_state) // TODO: Implement
	{
		std::cout << "Method SpinSpace::RemoveSubspace(const spin_ptr&, const state_ptr&) is not implemented!" << std::endl;
		return false;
	}

	// Remove the given spin and all spins that are coupled to it through the given interaction
	bool SpinSpace::RemoveSubspace(const spin_ptr &_spin, const interaction_ptr &_interaction) // TODO: Implement
	{
		std::cout << "Method SpinSpace::RemoveSubspace(const spin_ptr&, const interaction_ptr&) is not implemented!" << std::endl;
		return false;
	}

	// Removes all spins
	void SpinSpace::ClearSpins()
	{
		this->spins.clear();
	}
	// -----------------------------------------------------
	// Spin Management: Contains-methods
	// -----------------------------------------------------
	// Checks if the spin is contained within the spin space
	bool SpinSpace::Contains(const spin_ptr &_spin) const
	{
		if (std::find(this->spins.cbegin(), this->spins.cend(), _spin) != this->spins.cend())
			return true;

		return false;
	}

	// Checks whether all of the spins are contained within the spin space
	bool SpinSpace::Contains(const std::vector<spin_ptr> &_list) const // TODO: Implement more efficient version of the method
	{
		bool found_all = true;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			found_all &= this->Contains(*i);
		return found_all;
	}

	// Checks whether all of the spins in the given spin system are contained in the spin space
	bool SpinSpace::Contains(const std::shared_ptr<SpinSystem> &_system) const
	{
		if (_system == nullptr)
			return false;

		return this->Contains(_system->Spins());
	}

	// Same as above
	bool SpinSpace::Contains(const SpinSystem &_system) const
	{
		return this->Contains(_system.Spins());
	}

	bool SpinSpace::Contains(const std::string SpinObject) const
	{
		for (auto i = spins.cbegin(); i != spins.cend(); i++)
		{
			if(i->get()->Properties()->Name() == SpinObject)
			{
				return true;
			}
		}
		return false;
	}

	// Checks whether the spins space contains the given spin and all spins entangled to it in the given state
	bool SpinSpace::ContainsSubspace(const spin_ptr &_spin, const state_ptr &_state) const
	{
		std::cout << "Method SpinSpace::ContainsSubspace(const spin_ptr&, const state_ptr&) is not implemented!" << std::endl;
		return false;
	}

	// Checks whether the spin space contains the given spin and any spins coupled to it through the given interaction
	bool SpinSpace::ContainsSubspace(const spin_ptr &_spin, const interaction_ptr &_interaction) const
	{
		std::cout << "Method SpinSpace::ContainsSubspace(const spin_ptr&, const state_ptr&) is not implemented!" << std::endl;
		return false;
	}
	// -----------------------------------------------------
	// Interaction Management
	// -----------------------------------------------------
	// Add interaction if not already contained by the SpinSpace, and returns true if is was added
	bool SpinSpace::Add(const interaction_ptr &_interaction)
	{
		if (_interaction == nullptr)
			return false;

		if (std::find(this->interactions.cbegin(), this->interactions.cend(), _interaction) != this->interactions.cend())
			return false;

		this->interactions.push_back(_interaction);
		return true;
	}

	// Adds all the interactions from a given list, without allowing any duplicates. Returns true if any elements were added
	// NOTE: This method may change the order of the interactions in the vector if and only if new elements were inserted. TODO: Fix this before it causes confusion.
	bool SpinSpace::Add(const std::vector<interaction_ptr> &_list)
	{
		// Merge the two lists without keeping any duplicates
		std::vector<interaction_ptr> tmpvec;												 // Make a temporary list...
		tmpvec.reserve(this->interactions.size() + _list.size());							 // ...and reserve memory to hold the following three inserts
		tmpvec.insert(tmpvec.end(), _list.cbegin(), _list.cend());							 // Insert all new elements from the list
		tmpvec.insert(tmpvec.end(), this->interactions.cbegin(), this->interactions.cend()); // Insert all previous elements
		std::sort(tmpvec.begin(), tmpvec.end());											 // Sort the elements
		auto newEnd = std::unique(tmpvec.begin(), tmpvec.end());							 // Remove duplicates (requires list to be sorted)
		tmpvec.resize(std::distance(tmpvec.begin(), newEnd));								 // Shrink vector size to get rid of removed objects (std::unique cannot do that)

		// Check if any interactions were added (if the new collection contains more elements than the old)
		if (tmpvec.size() > this->interactions.size())
		{
			this->interactions = tmpvec;
			return true;
		}

		return false;
	}

	// Remove a single interaction
	bool SpinSpace::Remove(const interaction_ptr &_interaction)
	{
		auto i = std::find(this->interactions.cbegin(), this->interactions.cend(), _interaction);
		if (i != this->interactions.cend())
		{
			this->interactions.erase(i);
			return true;
		}

		return false;
	}

	// Remove a list of interactions
	bool SpinSpace::Remove(const std::vector<interaction_ptr> &_list) // TODO: Implement more efficient version of the method
	{
		bool removed_any = false;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			removed_any |= this->Remove(*i);
		return removed_any;
	}

	// Removes all interactions
	void SpinSpace::ClearInteractions()
	{
		this->interactions.clear();
	}

	// Checks whether the interaction is included in the spin space
	bool SpinSpace::Contains(const interaction_ptr &_interaction) const
	{
		if (std::find(this->interactions.cbegin(), this->interactions.cend(), _interaction) != this->interactions.cend())
			return true;

		return false;
	}

	// Checks whether every interaction in the list is contained in the spin space
	bool SpinSpace::Contains(const std::vector<interaction_ptr> &_list) const // TODO: Implement more efficient version of the method
	{
		bool found_all = true;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			found_all &= this->Contains(*i);
		return found_all;
	}

	// -----------------------------------------------------
	// Pulse Management
	// -----------------------------------------------------
	// Add pulse if not already contained by the SpinSpace, and returns true if is was added
	bool SpinSpace::Add(const pulse_ptr &_pulse)
	{
		if (_pulse == nullptr)
			return false;

		if (std::find(this->pulses.cbegin(), this->pulses.cend(), _pulse) != this->pulses.cend())
			return false;

		this->pulses.push_back(_pulse);
		return true;
	}

	// Adds all the pulses from a given list, without allowing any duplicates. Returns true if any elements were added
	// NOTE: This method may change the order of the pulses in the vector if and only if new elements were inserted. TODO: Fix this before it causes confusion.
	bool SpinSpace::Add(const std::vector<pulse_ptr> &_list)
	{
		// Merge the two lists without keeping any duplicates
		std::vector<pulse_ptr> tmpvec;												 // Make a temporary list...
		tmpvec.reserve(this->pulses.size() + _list.size());							 // ...and reserve memory to hold the following three inserts
		tmpvec.insert(tmpvec.end(), _list.cbegin(), _list.cend());							 // Insert all new elements from the list
		tmpvec.insert(tmpvec.end(), this->pulses.cbegin(), this->pulses.cend()); // Insert all previous elements
		std::sort(tmpvec.begin(), tmpvec.end());											 // Sort the elements
		auto newEnd = std::unique(tmpvec.begin(), tmpvec.end());							 // Remove duplicates (requires list to be sorted)
		tmpvec.resize(std::distance(tmpvec.begin(), newEnd));								 // Shrink vector size to get rid of removed objects (std::unique cannot do that)

		// Check if any pulses were added (if the new collection contains more elements than the old)
		if (tmpvec.size() > this->pulses.size())
		{
			this->pulses = tmpvec;
			return true;
		}

		return false;
	}

	// Remove a single pulse
	bool SpinSpace::Remove(const pulse_ptr &_pulse)
	{
		auto i = std::find(this->pulses.cbegin(), this->pulses.cend(), _pulse);
		if (i != this->pulses.cend())
		{
			this->pulses.erase(i);
			return true;
		}

		return false;
	}

	// Remove a list of pulses
	bool SpinSpace::Remove(const std::vector<pulse_ptr> &_list) // TODO: Implement more efficient version of the method
	{
		bool removed_any = false;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			removed_any |= this->Remove(*i);
		return removed_any;
	}

	// Removes all pulses
	void SpinSpace::ClearPulses()
	{
		this->pulses.clear();
	}

	// Checks whether the pulse is included in the spin space
	bool SpinSpace::Contains(const pulse_ptr &_pulse) const
	{
		if (std::find(this->pulses.cbegin(), this->pulses.cend(), _pulse) != this->pulses.cend())
			return true;

		return false;
	}

	// Checks whether every pulse in the list is contained in the spin space
	bool SpinSpace::Contains(const std::vector<pulse_ptr> &_list) const // TODO: Implement more efficient version of the method
	{
		bool found_all = true;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			found_all &= this->Contains(*i);
		return found_all;
	}

	// -----------------------------------------------------
	// Transition Management
	// -----------------------------------------------------
	// Add transition if not already contained by the SpinSpace, and returns true if is was added
	bool SpinSpace::Add(const transition_ptr &_transition)
	{
		if (_transition == nullptr)
			return false;

		if (std::find(this->transitions.cbegin(), this->transitions.cend(), _transition) != this->transitions.cend())
			return false;

		this->transitions.push_back(_transition);
		return true;
	}

	// Adds all the transitions from a given list, without allowing any duplicates. Returns true if any elements were added
	// NOTE: This method may change the order of the interactions in the vector if and only if new elements were inserted. TODO: Fix this before it causes confusion.
	bool SpinSpace::Add(const std::vector<transition_ptr> &_list)
	{
		// Merge the two lists without keeping any duplicates
		std::vector<transition_ptr> tmpvec;												   // Make a temporary list...
		tmpvec.reserve(this->transitions.size() + _list.size());						   // ...and reserve memory to hold the following three inserts
		tmpvec.insert(tmpvec.end(), _list.cbegin(), _list.cend());						   // Insert all new elements from the list
		tmpvec.insert(tmpvec.end(), this->transitions.cbegin(), this->transitions.cend()); // Insert all previous elements
		std::sort(tmpvec.begin(), tmpvec.end());										   // Sort the elements
		auto newEnd = std::unique(tmpvec.begin(), tmpvec.end());						   // Remove duplicates (requires list to be sorted)
		tmpvec.resize(std::distance(tmpvec.begin(), newEnd));							   // Shrink vector size to get rid of removed objects (std::unique cannot do that)

		// Check if any transitions were added (if the new collection contains more elements than the old)
		if (tmpvec.size() > this->transitions.size())
		{
			this->transitions = tmpvec;
			return true;
		}

		return false;
	}

	// Remove a single transition
	bool SpinSpace::Remove(const transition_ptr &_transition)
	{
		auto i = std::find(this->transitions.cbegin(), this->transitions.cend(), _transition);
		if (i != this->transitions.cend())
		{
			this->transitions.erase(i);
			return true;
		}

		return false;
	}

	// Remove a list of transitions
	bool SpinSpace::Remove(const std::vector<transition_ptr> &_list) // TODO: Implement more efficient version of the method
	{
		bool removed_any = false;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			removed_any |= this->Remove(*i);
		return removed_any;
	}

	// Removes all transitions
	void SpinSpace::ClearTransitions()
	{
		this->transitions.clear();
	}

	// Checks whether the transition is included in the spin space
	bool SpinSpace::Contains(const transition_ptr &_transition) const
	{
		if (std::find(this->transitions.cbegin(), this->transitions.cend(), _transition) != this->transitions.cend())
			return true;

		return false;
	}

	// Checks whether every transition in the list is contained in the spin space
	bool SpinSpace::Contains(const std::vector<transition_ptr> &_list) const // TODO: Implement more efficient version of the method
	{
		bool found_all = true;
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
			found_all &= this->Contains(*i);
		return found_all;
	}
}
