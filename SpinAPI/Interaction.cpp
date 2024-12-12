/////////////////////////////////////////////////////////////////////////
// Interaction class (SpinAPI Module)
// ------------------
// Base class for interactions.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <memory>
#include <iostream>
#include <Spin.h>
#include "ObjectParser.h"
#include "Interaction.h"

namespace SpinAPI
{
	// -----------------------------------------------------
	// Interaction Constructors and Destructor
	// -----------------------------------------------------
	// The constructor sets up the interaction parameters, but
	// the spin groups are read in the method ParseSpinGroups instead.
	Interaction::Interaction(std::string _name, std::string _contents) : properties(std::make_shared<MSDParser::ObjectParser>(_name, _contents)), couplingTensor(nullptr), tdTensor(3,3, arma::fill::zeros),
																		 field({0, 0, 0}), dvalue(0.0), evalue(0.0), group1(), group2(), type(InteractionType::Undefined), fieldType(InteractionFieldType::Static), prefactor(1.0), addCommonPrefactor(true), ignoreTensors(false), isValid(true),
																		 trjHasTime(false), trjHasField(false), trjHasPrefactor(false), trjHasTensor(false), trjTime(0), trjFieldX(0), trjFieldY(0), trjFieldZ(0), trjPrefactor(0),
																		 tdFrequency(1.0), tdPhase(0.0), tdAxis("0 0 1"), tdPerpendicularOscillation(false), tdInitialField({0, 0, 0}), tdInitialTensor(3,3, arma::fill::zeros)
	{//trjMatXX(0),trjMatXY(0),trjMatXZ(0),trjMatYX(0),trjMatYY(0), trjMatYZ(0),trjMatZX(0),trjMatZY(0),trjMatZZ(0),
		// Is a trajectory specified?
		std::string str;
		if (this->properties->Get("trajectory", str))
		{
			// Get the directory
			std::string directory;
			if (!this->properties->Get("FileReaderDirectory", directory))
			{
				std::cout << "Warning: Failed to obtain input file directory while initializing Interaction " << this->Name() << "!\n";
				std::cout << "         This may lead to problems when trying to localize the trajectory file." << std::endl;
			}

			// Load the trajectory, and check whether something was loaded
			if (this->trajectory.Load(directory + str))
			{
				if (this->trajectory.Length() > 0)
				{
					// Get column numbers
					this->trjHasTime = this->trajectory.HasColumn("time", trjTime);
					this->trjHasPrefactor = this->trajectory.HasColumn("prefactor", trjPrefactor);
					this->trjHasField = this->trajectory.HasColumn("field.x", trjFieldX);
					this->trjHasField &= this->trajectory.HasColumn("field.y", trjFieldY);
					this->trjHasField &= this->trajectory.HasColumn("field.z", trjFieldZ);

					// this->trjHasTensor = this->trajectory.HasColumn("mat.xx", trjMatXX);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.xy", trjMatXY);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.xz", trjMatXZ);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.yx", trjMatYX);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.yy", trjMatYY);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.yz", trjMatYZ);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.zx", trjMatZX);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.zy", trjMatZY);
					// this->trjHasTensor &= this->trajectory.HasColumn("mat.zz", trjMatZZ);
				}
			}
			else
			{
				std::cout << "ERROR: Failed to load trajectory \"" << (directory + str) << "\" for Interaction object \"" << this->Name() << "\"!" << std::endl;
			}
		}

		// Get the type of the interaction
		if (this->properties->Get("type", str))
		{
			if (str.compare("onespin") == 0 || str.compare("singlespin") == 0 || str.compare("zeeman") == 0)
			{
				// A onespin interaction should have a field attached, either directly or through a trajectory
				arma::vec inField = arma::zeros<arma::vec>(3);
				if (this->properties->Get("field", inField) || this->trjHasField)
				{
					this->field = inField;
					this->type = InteractionType::SingleSpin;
				}
			}
			else if (str.compare("twospin") == 0 || str.compare("doublespin") == 0 || str.compare("hyperfine") == 0 || str.compare("dipole") == 0)
			{
				this->type = InteractionType::DoubleSpin;
			}
			else if (str.compare("exchange") == 0)
			{
				this->type = InteractionType::Exchange;
			}
			else if (str.compare("zfs") == 0)
			{
				this->type = InteractionType::Zfs;

				double indvalue, inevalue;
				this->Properties()->Get("dvalue", indvalue);
				this->Properties()->Get("evalue", inevalue);

				this->dvalue = indvalue;
				this->evalue = inevalue;
			}
		}

		// If we have a valid interaction type, read the other parameters
		if (this->type != InteractionType::Undefined)
		{
			Tensor inTensor(0);
			if (this->properties->Get("tensor", inTensor))
			{
				this->couplingTensor = std::make_shared<Tensor>(inTensor);
			}

			this->properties->Get("prefactor", this->prefactor);
			this->properties->Get("commonprefactor", this->addCommonPrefactor);
			this->properties->Get("ignoretensors", this->ignoreTensors);
		}

		// One-spin interactions can have a time-dependent field
		if (this->type == InteractionType::SingleSpin)
		{
			// Do we have a trajectory entry for field? Note that the trajectory may not have a time column, in which case there will be no time dependence
			if (this->trjHasField)
			{
				this->fieldType = InteractionFieldType::Trajectory;
			}

			if (this->properties->Get("fieldtype", str) || this->properties->Get("timedependence", str))
			{
				if (this->trjHasField)
				{
					std::cout << "Warning: Ignored fieldtype for Interaction \"" << this->Name() << "\"! Using trajectory instead." << std::endl;
				}
				// Figure out what type of time-dependence should be used
				else if (str.compare("linearpolarized") == 0 || str.compare("linearpolarization") == 0 || str.compare("oscillating") == 0)
				{
					this->fieldType = InteractionFieldType::LinearPolarization;
					this->properties->Get("frequency", this->tdFrequency);
					this->properties->Get("phase", this->tdPhase);
				}
				else if (str.compare("circularpolarized") == 0 || str.compare("circularpolarization") == 0)
				{
					this->fieldType = InteractionFieldType::CircularPolarization;
					this->properties->Get("frequency", this->tdFrequency);
					this->properties->Get("phase", this->tdPhase);
					this->properties->Get("axis", this->tdAxis);
					this->properties->Get("perpendicularoscillations", this->tdPerpendicularOscillation);
				}
				else
				{
					std::cout << "Warning: Unknown fieldtype for Interaction \"" << this->Name() << "\"! Assuming static field." << std::endl;
				}
			}

			// Save the initial field - used only for time-dependent interactions where Actions cannot change the field
			this->tdInitialField = this->field;

			// Set the time to 0 by default
			if (this->HasFieldTimeDependence())
				this->SetTime(0.0);
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////PIP ADDITIONS /////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// Double-spin interactions can have a time-dependent tensor
		if (this->type == InteractionType::DoubleSpin)
		{

			// Do we have a trajectory entry for field? Note that the trajectory may not have a time column, in which case there will be no time dependence
			if (this->trjHasTensor)
			{
				this->tensorType = InteractionTensorType::Trajectory;
			}

			if (this->properties->Get("fieldtype", str) || this->properties->Get("timedependence", str))
			{
				if (this->trjHasTensor)
				{
					std::cout << "Warning: Ignored fieldtype for Interaction \"" << this->Name() << "\"! Using trajectory instead." << std::endl;
				}
				// Figure out what type of time-dependence should be used
				else if (str.compare("SinMat") == 0)
				{
					this->tensorType = InteractionTensorType::SinMat;
					this->properties->Get("frequency", this->tdFrequency);
					this->properties->Get("phase", this->tdPhase);
				}
				else
				{
					std::cout << "Warning: Unknown fieldtype for Interaction \"" << this->Name() << "\"! Assuming static field." << std::endl;
				}
			
			}

			// Save the initial tensor - used only for time-dependent interactions where Actions cannot change the field
			arma::mat zeros(3,3, arma::fill::zeros);
			this->tdInitialTensor = zeros;//this->couplingTensor;

			//Set the time to 0 by default
			if (this->HasFieldTimeDependence())
				this->SetTime(0.0);

		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////PIP ADDITIONS /////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	Interaction::Interaction(const Interaction &_interaction) : properties(_interaction.properties), couplingTensor(_interaction.couplingTensor), field(_interaction.field), dvalue(_interaction.dvalue), evalue(_interaction.evalue),
																group1(_interaction.group1), group2(_interaction.group2), type(_interaction.type), fieldType(_interaction.fieldType),
																prefactor(_interaction.prefactor), addCommonPrefactor(_interaction.addCommonPrefactor), ignoreTensors(_interaction.ignoreTensors), isValid(_interaction.isValid),
																trjHasTime(_interaction.trjHasTime), trjHasField(_interaction.trjHasField), trjHasPrefactor(_interaction.trjHasPrefactor), trjHasTensor(_interaction.trjHasTensor),
																trjTime(_interaction.trjTime), trjFieldX(_interaction.trjFieldX), trjFieldY(_interaction.trjFieldY), trjFieldZ(_interaction.trjFieldZ),
																trjPrefactor(_interaction.trjPrefactor), tdFrequency(_interaction.tdFrequency), tdPhase(_interaction.tdPhase), tdAxis(_interaction.tdAxis),
																tdPerpendicularOscillation(_interaction.tdPerpendicularOscillation), tdInitialField(_interaction.tdInitialField), tdInitialTensor(_interaction.tdInitialTensor)
	{
	}

	Interaction::~Interaction()
	{
	}
	// -----------------------------------------------------
	// Operators
	// -----------------------------------------------------
	const Interaction &Interaction::operator=(const Interaction &_interaction)
	{
		this->properties = std::make_shared<MSDParser::ObjectParser>(*(_interaction.properties));
		this->couplingTensor = _interaction.couplingTensor;
		this->field = _interaction.field;
		this->dvalue = _interaction.dvalue;
		this->evalue = _interaction.evalue;
		this->type = _interaction.type;
		this->fieldType = _interaction.fieldType;
		this->prefactor = _interaction.prefactor;
		this->addCommonPrefactor = _interaction.addCommonPrefactor;
		this->ignoreTensors = _interaction.ignoreTensors;
		this->isValid = _interaction.isValid;

		this->trjHasTime = _interaction.trjHasTime;
		this->trjHasField = _interaction.trjHasField;
		this->trjHasTensor = _interaction.trjHasTensor;
		this->trjHasPrefactor = _interaction.trjHasPrefactor;
		this->trjTime = _interaction.trjTime;
		this->trjFieldX = _interaction.trjFieldX;
		this->trjFieldY = _interaction.trjFieldY;
		this->trjFieldZ = _interaction.trjFieldZ;
		this->trjPrefactor = _interaction.trjPrefactor;

		this->tdFrequency = _interaction.tdFrequency;
		this->tdPhase = _interaction.tdPhase;
		this->tdAxis = _interaction.tdAxis;
		this->tdPerpendicularOscillation = _interaction.tdPerpendicularOscillation;
		this->tdInitialField = _interaction.tdInitialField;
		this->tdInitialTensor = _interaction.tdInitialTensor;
		this->tdTensor = _interaction.tdTensor;

		return (*this);
	}
	// -----------------------------------------------------
	// Public methods
	// -----------------------------------------------------
	std::string Interaction::Name() const
	{
		return this->properties->Name();
	}

	bool Interaction::IsValid() const
	{
		if(!isValid) { return false; } //by default this function should never return at this point; the only time it should ever return is if the interaction hasn't been assigned a subsystem in a task where subsystems are used 
		
		if (this->type == InteractionType::SingleSpin && !this->group1.empty())
			return true;
		else if (this->type == InteractionType::DoubleSpin && !this->group1.empty() && !this->group2.empty())
			return true;
		else if (this->type == InteractionType::Exchange && !this->group1.empty() && !this->group2.empty())
			return true;
		else if (this->type == InteractionType::Zfs && !this->group1.empty())
			return true;

		return false;
	}
	// -----------------------------------------------------
	// Property methods
	// -----------------------------------------------------
	// Returns the field vector
	const arma::vec Interaction::Field() const
	{
		return this->field;
	}

	// Returns the D value for Zfs
	const double Interaction::Dvalue() const
	{
		return this->dvalue;
	}

	// Returns the E value for Zfs
	const double Interaction::Evalue() const
	{
		return this->evalue;
	}

	// Returns the prefactor value
	const double Interaction::Prefactor() const
	{
		return this->prefactor;
	}

	// Checks whether the interaction is time-dependent
	bool Interaction::HasFieldTimeDependence() const
	{
		if (this->fieldType != InteractionFieldType::Static)
			return true;

		return false;
	}

	// Checks whether the interaction is time-dependent
	bool Interaction::HasTimeDependence() const
	{
		if ((this->type == InteractionType::SingleSpin && this->HasFieldTimeDependence()) || (this->trjHasTime && this->trjHasPrefactor))
			return true;

		if (this->trjHasTime)
		{
			return true;
		}

		return false;
	}

	//TODO : ADD IN TENSOR TIMEDEPENDENCE CONDITIONAL METHODS
	

	// -----------------------------------------------------
	// Trajectory-related methods
	// -----------------------------------------------------
	// Returns the length of the trajectory, if any
	unsigned int Interaction::TrajectoryLength() const
	{
		return this->trajectory.Length();
	}

	// Set the current state of the Interaction to a specific step in the trajectory (row in the trajectory file)
	bool Interaction::SetTrajectoryStep(unsigned int _step)
	{
		// Make sure that we have a trajectory
		if (this->trajectory.Length() < 1 || _step >= this->trajectory.Length())
		{
			if (this->couplingTensor != nullptr)
				return this->couplingTensor->SetTrajectoryStep(_step); // The tensor might also have a trajectory, so don't return false without asking the tensor...
			else
				return false;
		}

		// Check whether the requested trajectory step is valid
		if (_step < this->trajectory.Length())
		{
			// Get the prefactor at the requested trajectory step
			if (this->trjHasPrefactor)
				this->prefactor = this->trajectory.Get(_step, this->trjPrefactor);

			// Get the field at the requested trajectory step
			if (this->trjHasField && this->fieldType == InteractionFieldType::Trajectory)
			{
				arma::vec tmp = arma::zeros<arma::vec>(3);
				tmp(0) = this->trajectory.Get(_step, this->trjFieldX);
				tmp(1) = this->trajectory.Get(_step, this->trjFieldY);
				tmp(2) = this->trajectory.Get(_step, this->trjFieldZ);
				this->field = tmp;
			}

			// Also set the trajectory step for the tensor, if possible
			// TODO: Consider what should happen if tensor trajectory length is different from interaction object trajectory length
			if (this->couplingTensor != nullptr)
				this->couplingTensor->SetTrajectoryStep(_step);
		}

		return true;
	}

	// Set the current state of the tensor to a specific time defined either in the trajectory (requires a "time" column in the trajectory) or the time-dependence function of the field
	bool Interaction::SetTime(double _time)
	{
		// Set the time for the tensor; the tensor trajectory is independent of the trajectory/time dependence of the interaction object
		if (this->couplingTensor != nullptr)
			this->couplingTensor->SetTime(_time);

		// If we have a trajectory, which has a "time" column
		if (this->trajectory.Length() > 0 && this->trjHasTime)
		{
			// Get the row corresponding to the time right after (or at) the requested time
			unsigned int row = 0;
			double timestepAbove = 0.0;
			if (!this->trajectory.FirstRowEqGreaterThan(_time, this->trjTime, row, timestepAbove))
			{
				// If no such time was found, i.e. we are past the last element of the trajectory, use the last step of the trajectory
				this->SetTrajectoryStep(this->trajectory.Length() - 1);
			}
			else if (row == 0)
			{
				// If the specified time is before the first element of the trajectory, use the first step of the trajectory
				this->SetTrajectoryStep((unsigned int)0);
			}
			else
			{
				// Prepare a linear interpolation
				double timestepBelow = this->trajectory.Get(row - 1, trjTime);			  // Get time at previous step
				double l = 1 - (timestepAbove - _time) / (timestepAbove - timestepBelow); // Get linear interpolation factor, i.e. a number between 0 and 1, where 0 = previous time step and 1 = next time step in the trajectory

				// Get a linear interpolation for the prefactor
				if (this->trjHasPrefactor)
					this->prefactor = this->trajectory.Get(row, this->trjPrefactor) * l + this->trajectory.Get(row - 1, this->trjPrefactor) * (1 - l);

				// Get a linear interpolation for the field vector
				if (this->trjHasField && this->fieldType == InteractionFieldType::Trajectory)
				{
					arma::vec tmp = arma::zeros<arma::vec>(3);
					tmp(0) = this->trajectory.Get(row, this->trjFieldX) * l + this->trajectory.Get(row - 1, this->trjFieldX) * (1 - l);
					tmp(1) = this->trajectory.Get(row, this->trjFieldY) * l + this->trajectory.Get(row - 1, this->trjFieldY) * (1 - l);
					tmp(2) = this->trajectory.Get(row, this->trjFieldZ) * l + this->trajectory.Get(row - 1, this->trjFieldZ) * (1 - l);
					this->field = tmp;
				}
			}
		}

		// Apply the right time-dependence function if not using a trajectory or static field
		if (this->fieldType == InteractionFieldType::LinearPolarization)
			this->field = FieldTimeDependenceLinearPolarization(this->tdInitialField, _time, this->tdFrequency, this->tdPhase);
		else if (this->fieldType == InteractionFieldType::CircularPolarization)
			this->field = FieldTimeDependenceCircularPolarization(this->tdInitialField, _time, this->tdFrequency, this->tdPhase, this->tdAxis, this->tdPerpendicularOscillation);
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////// TENSOR TIMEDEP FUNCTIONS IN HERE //////////////////////////////////////
		else if (this->tensorType == InteractionTensorType::SinMat){
			arma::mat td_tensor = TensorTimeDependenceSinMat(this->tdInitialTensor, _time, this->tdFrequency, this->tdPhase);
			// Tensor cp_tensor = Tensor(td_tensor);
			this->couplingTensor = std::make_shared<Tensor>(td_tensor);
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		return true;
	}

	// -----------------------------------------------------
	// Access to custom properties
	// -----------------------------------------------------
	std::shared_ptr<const MSDParser::ObjectParser> Interaction::Properties() const
	{
		return this->properties;
	}

	// -----------------------------------------------------
	// Public method to read spins into group1 and group2,
	// returns false if some spins were not found
	// -----------------------------------------------------
	bool Interaction::ParseSpinGroups(const std::vector<spin_ptr> &_spinlist)
	{
		bool createdSpinLists = false;

		if (this->type == InteractionType::SingleSpin || this->type == InteractionType::Zfs)
		{
			// Attempt to get a list of spins from the input file
			std::string str;
			if (!this->properties->Get("spins", str) && !this->properties->Get("group1", str))
				return false;

			createdSpinLists = this->AddSpinList(str, _spinlist, this->group1);
		}
		else if (this->type == InteractionType::DoubleSpin || this->type == InteractionType::Exchange)
		{
			// Attempt to get a list of spins from the input file
			std::string str1;
			std::string str2;
			if (!this->properties->Get("group1", str1) || !this->properties->Get("group2", str2))
				return false;

			createdSpinLists = this->AddSpinList(str1, _spinlist, this->group1);
			createdSpinLists &= this->AddSpinList(str2, _spinlist, this->group2, &(this->group1)); // Cross-check with group1 as a spin can be in only one of the groups
		}

		// Check whether we were successful
		if (!createdSpinLists)
		{
			this->group1.clear();
			this->group2.clear();
			return false;
		}

		return true;
	}

	// Splits the string into spin names, find the corresponding spins in _spinlist and adds them to _group
	// Private helper method for the ParseSpinGroups method
	bool Interaction::AddSpinList(const std::string &_names, const std::vector<spin_ptr> &_spinlist, std::vector<spin_ptr> &_group, const std::vector<spin_ptr> *_crossCheck)
	{
		// Iterate through the spin name (while-loop) with a ',' delimiter
		std::stringstream ss(_names);
		bool spin_was_found;
		std::string str;
		while (std::getline(ss, str, ','))
		{
			spin_was_found = false;

			// Get the spin with the name held by str
			for (auto i = _spinlist.cbegin(); i != _spinlist.cend(); i++)
			{
				if ((*i)->Name().compare(str) == 0)
				{
					// Perform crosscheck if requested
					if (_crossCheck != nullptr && std::find((*_crossCheck).cbegin(), (*_crossCheck).cend(), (*i)) != (*_crossCheck).cend())
						return false;

					// Make sure we don't have duplicates in the spin group
					if (std::find(_group.cbegin(), _group.cend(), (*i)) != _group.cend())
						return false;

					// Add it to the group
					_group.push_back(*i);
					spin_was_found = true;
					break;
				}
			}

			if (!spin_was_found)
				return false;
		}

		return true;
	}
	// -----------------------------------------------------
	// Public spin list methods
	// -----------------------------------------------------
	// Returns a list of spins such that no spins in the list interacts with a spin outside the list (via the current Interaction object)
	std::vector<spin_ptr> Interaction::CompleteSet(const spin_ptr &_spin) const
	{
		std::vector<spin_ptr> result(1); // Only allocate space for a single spin_ptr by default

		// Only the DoubleSpin interaction type couples spins.
		// If the spin is found in either group, return a list containing the set-theoretic union of the groups
		if (this->type == InteractionType::DoubleSpin &&
			(std::find(this->group1.cbegin(), this->group1.cend(), _spin) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), _spin) != this->group2.cend()))
		{
			result.reserve(this->group1.size() + this->group2.size());				   // Reserve space for both groups to avoid more than 1 reallocation
			result.insert(result.begin(), this->group1.cbegin(), this->group1.cend()); // Insert at the beginning of the vector
			result.insert(result.end(), this->group2.cbegin(), this->group2.cend());   // Insert after the previously inserted spins
		}
		else if (this->type == InteractionType::Exchange &&
				 (std::find(this->group1.cbegin(), this->group1.cend(), _spin) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), _spin) != this->group2.cend()))
		{
			result.reserve(this->group1.size() + this->group2.size());				   // Reserve space for both groups to avoid more than 1 reallocation
			result.insert(result.begin(), this->group1.cbegin(), this->group1.cend()); // Insert at the beginning of the vector
			result.insert(result.end(), this->group2.cbegin(), this->group2.cend());   // Insert after the previously inserted spins
		}
		else if (this->type == InteractionType::Zfs &&
				 (std::find(this->group1.cbegin(), this->group1.cend(), _spin) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), _spin) != this->group2.cend()))
		{
			result.reserve(this->group1.size() + this->group2.size());				   // Reserve space for both groups to avoid more than 1 reallocation
			result.insert(result.begin(), this->group1.cbegin(), this->group1.cend()); // Insert at the beginning of the vector
			result.insert(result.end(), this->group2.cbegin(), this->group2.cend());   // Insert after the previously inserted spins
		}
		else
		{
			// If the spin is not in either group, the spin itself forms a complete set with respect to the current interaction object
			result.push_back(_spin);
		}

		return result;
	}

	// Extends the list of spins such that no spins in the list interacts with a spin outside the list (via the current Interaction object)
	// NOTE: The order of the spins in the list is changed if and only if the list is extended.
	bool Interaction::CompleteSet(std::vector<spin_ptr> &_list) const
	{
		bool was_extended = false;

		// Only the DoubleSpin interaction type couples spins
		if (this->type == InteractionType::DoubleSpin || this->type == InteractionType::Exchange)
		{
			bool containsInteractingSpins = false;

			// Loop through all elements in the list
			for (auto i = _list.cbegin(); i != _list.cend(); i++)
			{
				// Check whether any of the spins in the list are coupled to another spin via the current Interaction object
				if (std::find(this->group1.cbegin(), this->group1.cend(), (*i)) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), (*i)) != this->group2.cend())
				{
					containsInteractingSpins = true;
					break;
				}
			}

			// There are interacting spins in the list, so we may need to extend it
			if (containsInteractingSpins)
			{
				// First get a temporary list that is the set-theoretic completion of the given list
				std::vector<spin_ptr> tmpvec;											  // Make a temporary list...
				tmpvec.reserve(_list.size() + this->group1.size() + this->group2.size()); // ...and reserve memory to hold the following three inserts
				tmpvec.insert(tmpvec.end(), _list.cbegin(), _list.cend());				  // Insert all elements from the list
				tmpvec.insert(tmpvec.end(), this->group1.cbegin(), this->group1.cend());  // Insert all elements from group1
				tmpvec.insert(tmpvec.end(), this->group2.cbegin(), this->group2.cend());  // Insert all elements from group2
				std::sort(tmpvec.begin(), tmpvec.end());								  // Sort the elements
				auto newEnd = std::unique(tmpvec.begin(), tmpvec.end());				  // Remove duplicates (requires list to be sorted)
				tmpvec.resize(std::distance(tmpvec.begin(), newEnd));					  // Shrink vector size to get rid of removed objects (std::unique cannot do that)

				// Then check if we added new elements in the process
				if (tmpvec.size() > _list.size())
				{
					// Replace the list with the temporary one, if the temporary list was an extension
					_list = tmpvec;
					was_extended = true;
				}
			}
		}

		return was_extended;
	}

	// Public method to check whether a list of spins is complete,
	// i.e. does not interact with any spins outside the list
	bool Interaction::IsComplete(const std::vector<spin_ptr> &_list) const
	{
		unsigned int group1Elements = 0;
		unsigned int group2Elements = 0;

		// Loop through all elements in the list
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
		{
			// Count number of elements from group1 found in the list
			if (std::find(this->group1.cbegin(), this->group1.cend(), (*i)) != this->group1.cend())
				++group1Elements;

			// Count number of elements from group2 found in the list
			if (std::find(this->group2.cbegin(), this->group2.cend(), (*i)) != this->group2.cend())
				++group2Elements;
		}

		// If an element from group1 is present, all elements from group2 must be present, and vice versa
		if (group1Elements > 0 && group2Elements < this->group2.size())
			return false;
		else if (group2Elements > 0 && group1Elements < this->group1.size())
			return false;

		return true;
	}
	// -----------------------------------------------------
	// Methods to create ActionTarget objects
	// -----------------------------------------------------
	// Create ActionVectors
	std::vector<RunSection::NamedActionVector> Interaction::CreateActionVectors(const std::string &_system)
	{
		std::vector<RunSection::NamedActionVector> vectors;

		if (this->IsValid())
		{
			if (this->type == InteractionType::SingleSpin)
			{
				// The field vector
				RunSection::ActionVector fieldVector = RunSection::ActionVector(this->field, &CheckActionVectorInteractionField, this->HasFieldTimeDependence());
				vectors.push_back(RunSection::NamedActionVector(_system + "." + this->Name() + ".field", fieldVector));

				// Allows changing the base field used by time-dependence functions
				if (this->HasFieldTimeDependence() && this->fieldType != InteractionFieldType::Trajectory)
				{
					RunSection::ActionVector initialfieldVector = RunSection::ActionVector(this->tdInitialField, &CheckActionVectorInteractionField);
					vectors.push_back(RunSection::NamedActionVector(_system + "." + this->Name() + ".basefield", initialfieldVector));
				}

				// The axis of circular polarized oscillating fields
				if (this->fieldType == InteractionFieldType::CircularPolarization)
				{
					RunSection::ActionVector oscaxisVector = RunSection::ActionVector(this->tdAxis, &CheckActionVectorInteractionField);
					vectors.push_back(RunSection::NamedActionVector(_system + "." + this->Name() + ".axis", oscaxisVector));
				}
			}
		}

		return vectors;
	}

	// Create ActionScalars
	std::vector<RunSection::NamedActionScalar> Interaction::CreateActionScalars(const std::string &_system)
	{
		std::vector<RunSection::NamedActionScalar> scalars;

		if (this->IsValid())
		{
			// We should always have a scalar for the prefactor
			RunSection::ActionScalar prefactorScalar = RunSection::ActionScalar(this->prefactor, &CheckActionScalarInteractionPrefactor, this->trjHasPrefactor);
			scalars.push_back(RunSection::NamedActionScalar(_system + "." + this->Name() + ".prefactor", prefactorScalar));

			// There may be additional scalars for singlespin interactions
			if (this->type == InteractionType::SingleSpin)
			{
				// Frequency and phase
				if (this->fieldType == InteractionFieldType::LinearPolarization || this->fieldType == InteractionFieldType::CircularPolarization)
				{
					RunSection::ActionScalar frequencyScalar = RunSection::ActionScalar(this->tdFrequency, nullptr);
					scalars.push_back(RunSection::NamedActionScalar(_system + "." + this->Name() + ".frequency", frequencyScalar));

					RunSection::ActionScalar phaseScalar = RunSection::ActionScalar(this->tdPhase, nullptr);
					scalars.push_back(RunSection::NamedActionScalar(_system + "." + this->Name() + ".phase", phaseScalar));
				}
			}
		}

		return scalars;
	}

	// Method that calls the methods to generate ActionVectors and ActionScalars and inserts them into the given collections
	void Interaction::GetActionTargets(std::vector<RunSection::NamedActionScalar> &_scalars, std::vector<RunSection::NamedActionVector> &_vectors, const std::string &_system)
	{
		// Get ActionTargets from private methods
		auto scalars = this->CreateActionScalars(_system);
		auto vectors = this->CreateActionVectors(_system);

		// Insert them
		_scalars.insert(_scalars.end(), scalars.begin(), scalars.end());
		_vectors.insert(_vectors.end(), vectors.begin(), vectors.end());

		// If we have a tensor assigned, also add the ActionTargets of the Tensor
		if (this->couplingTensor != nullptr)
			this->couplingTensor->GetActionTargets(_scalars, _vectors, _system + "." + this->Name());
	}
	// -----------------------------------------------------
	// Non-member non-friend methods
	// -----------------------------------------------------
	bool IsValid(const Interaction &_interaction)
	{
		return _interaction.IsValid();
	}

	bool IsStatic(const Interaction &_interaction)
	{
		return !(HasTrajectory(_interaction) | _interaction.HasTimeDependence());
	}

	bool HasTensor(const Interaction &_interaction)
	{
		return (_interaction.CouplingTensor() != nullptr);
	}

	bool HasTrajectory(const Interaction &_interaction)
	{
		if (_interaction.TrajectoryLength() > 0)
			return true;

		return false;
	}

	InteractionType Type(const Interaction &_interaction)
	{
		return _interaction.Type();
	}
	// -----------------------------------------------------
	// Non-member non-friend time-dependent field functions
	// -----------------------------------------------------
	// Linear polarized oscillating magnetic field
	arma::vec FieldTimeDependenceLinearPolarization(const arma::vec &_v, double _time, double _frequency, double _phase)
	{
		return (_v * cos(_frequency * _time + _phase));
	}

	// Circularly polarized oscillating magnetic field
	arma::vec FieldTimeDependenceCircularPolarization(const arma::vec &_v, double _time, double _frequency, double _phase, const arma::vec &_axis, bool _perpendicularOscillations)
	{
		// Get the angle for the rotation matrix
		double angle = _time * _frequency + _phase;

		// Helper matrix
		arma::mat W(3, 3, arma::fill::zeros);
		W(0, 1) = -_axis(2);
		W(0, 2) = _axis(1);
		W(1, 2) = -_axis(0);
		W(1, 0) = _axis(2);
		W(2, 0) = -_axis(1);
		W(2, 1) = _axis(0);

		// The Rodrigues rotation matrix
		arma::mat R(3, 3, arma::fill::eye);
		R += sin(angle) * W + (2 * sin(angle / 2.0) * sin(angle / 2.0) * W * W);

		if (_perpendicularOscillations)
		{
			// Get the projection of the field vector in the plane perpendicular to the axis, and apply the rotation
			return R * (_v - dot(_axis, _v) * _axis);
		}

		// Just apply the rotation matrix to the field without projecting onto perpendicular plane
		return R * _v;
	}

	//make a matrix with sinusoidal modulation on one tensor component
	arma::mat TensorTimeDependenceSinMat(arma::mat _m, double _time, double _frequency, double _phase){
		
		// std::cout << "RECOGNISED" << std::endl;
		// _m.trjMatXX *= cos(_frequency * _time + _phase);
		_m(0,0) = 1 * cos(_frequency * _time + _phase);
		return _m;
		//REDEFINE IN TERMS OF MOLSPIN TENSOR

	}



	// -----------------------------------------------------
	// Non-member non-friend ActionTarget Check functions
	// -----------------------------------------------------
	// Checks that a vector set by an Action is valid
	bool CheckActionVectorInteractionField(const arma::vec &_v)
	{
		// A field vector must have 3 components
		if (_v.n_elem != 3)
			return false;

		// Make sure we don't have invalid values
		if (_v.has_nan() || _v.has_inf())
			return false;

		return true;
	}

	// Make sure that the prefactor has a valid value (not NaN or infinite)
	bool CheckActionScalarInteractionPrefactor(const double &_d)
	{
		return std::isfinite(_d);
	}
	// -----------------------------------------------------
}
