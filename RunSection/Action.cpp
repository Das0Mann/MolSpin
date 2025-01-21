/////////////////////////////////////////////////////////////////////////
// Action implementation (RunSection module)
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ObjectParser.h"
#include "Action.h"

namespace RunSection
{
	// -----------------------------------------------------
	// Action Constructors and Destructor
	// -----------------------------------------------------
	Action::Action(const MSDParser::ObjectParser &_properties, const std::map<std::string, ActionScalar> &_scalars, const std::map<std::string, ActionVector> &_vectors)
		: properties(std::make_shared<MSDParser::ObjectParser>(_properties)), scalars(_scalars), vectors(_vectors), isValid(false), value(1.0), first(0), last(0), period(1)
	{
		this->properties->Get("value", this->value);
		this->properties->Get("first", this->first);
		this->properties->Get("last", this->last);
		this->properties->Get("period", this->period);
		
		//std::string Parallelize = "";
		//this->properties->Get("parallelize", Parallelize);
		//if(Parallelize.compare("true") == 0)
		//{
		//	this->m_Parallelize = true;
		//}
		//else
		//{
		//	this->m_Parallelize = false;
		//}
		
		std::string loop;
		this->properties->Get("loop", loop);
		
		if(loop.compare("true") == 0)
		{
			this->m_loop = true;
		}
		else
		{
			this->m_loop = false;
		}
	}

	Action::~Action()
	{
	}
	// -----------------------------------------------------
	// Action protected methods
	// -----------------------------------------------------
	// Sets a reference to the requested ActionScalar
	bool Action::Scalar(std::string _name, ActionScalar **_scalar) const
	{
		auto i = this->scalars.find(_name);
		if (i != this->scalars.cend())
		{
			// Set the pointer if any
			if (_scalar != nullptr)
				(*_scalar) = const_cast<ActionScalar *>(&(i->second));

			return true;
		}

		return false;
	}

	// Sets a reference to the requested ActionVector
	bool Action::Vector(std::string _name, ActionVector **_vector) const
	{
		auto i = this->vectors.find(_name);
		if (i != this->vectors.cend())
		{
			// Set the pointer if any
			if (_vector != nullptr)
				(*_vector) = const_cast<ActionVector *>(&(i->second));

			return true;
		}

		return false;
	}
	// -----------------------------------------------------
	// Action public methods
	// -----------------------------------------------------
	// The specified value for the action
	double Action::Value() const
	{
		return this->value;
	}

	// The first step that should trigger the action
	unsigned int Action::First() const
	{
		return this->first;
	}

	// The last step that should trigger the action
	unsigned int Action::Last() const
	{
		return this->last;
	}

	// The periodicity for calls to the action
	unsigned int Action::Period() const
	{
		return this->period;
	}

	// Use the action
	void Action::Step(unsigned int _currentStep)
	{
		// Make sure that the action is valid
		if (!this->isValid)
			return;

		///bool TopLevel = true;
		///if(m_LoopLevel == -1 || m_steps == -1)
		///{
		///	TopLevel = true;
		///}

		// Only use the action if we are at a step between "first" and "last", and take steps with the given periodicity
		if(!m_loop)
		{
			if (_currentStep < this->first || (_currentStep > this->last && this->last != 0) || ((_currentStep - this->first - 1) % this->period) != 0)
				return;

			//std::cout << "updating" << std::endl;
		}

		if(m_loop)
		{
			if(_currentStep == this->last + 1 && this->last != 0)
			{
				this->Reset();
				//std::cout << "resetting" << std::endl; //current code assumes starting from 0, commented code assumes startitng from 1
				int gap = this->last - this->first; //this-last - this-first
				this->first = _currentStep; //_currentStep
				this->last = this->first + gap; //this->first + gap 
				return;
			}

			if(((_currentStep - this->first -1) % this->period) != 0)
			{
				return;
			}

			if(_currentStep < this->first)
			{
				return;
			}
		}

		// Take a step with the action, and notify the user if it failed
		if (!this->DoStep())
			std::cerr << "ERROR: Failed to perform step with action \"" << this->Name() << "\"!" << std::endl;
	}

	// Prepares the action and checks whether it is valid
	bool Action::Validate()
	{
		// No need to create an action object if it is never used
		if (this->first > this->last && this->last != 0)
			return false;

		// Let derived classes validate their input
		this->isValid = this->DoValidate();
		return this->isValid;
	}

	// Checks whether the action is valid
	bool Action::IsValid() const
	{
		return this->isValid;
	}

	// Returns the name of the action object
	std::string Action::Name()
	{
		if (this->properties == nullptr)
			return "undefined";

		return this->properties->Name();
	}
	// -----------------------------------------------------
}
