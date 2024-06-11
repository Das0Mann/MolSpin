/////////////////////////////////////////////////////////////////////////
// Interaction class (SpinAPI Module)
// ------------------
// Base class for interactions.
// 
// Molecular Spin Dynamics Software - developed by Luca.
// (c) 2024 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <memory>
#include <iostream>
#include <Spin.h>
#include "ObjectParser.h"
#include "Pulse.h"

namespace SpinAPI
{
	// -----------------------------------------------------
	// Interaction Constructors and Destructor
	// -----------------------------------------------------
	// The constructor sets up the pulse parameters, but
	// the spin groups are read in the method ParseSpinGroups instead.
	Pulse::Pulse(std::string _name, std::string _contents)	: properties(std::make_shared<MSDParser::ObjectParser>(_name,_contents)), type(PulseType::Unspecified), group(), rotationaxis({0, 0, 0}), angle(0.0), pulsetime(0.001), field({0, 0, 0}) 
	{
        // Filling required parameter

		// Get rotation axis
		if(!this->properties->Get("rotationaxis", rotationaxis))
		{
			std::cout << "Warning: Failed to obtain input for rotationaxis." << std::endl;
			std::cout << "Using vector: " << rotationaxis << std::endl;
		}

		// Get angle
		if(!this->properties->Get("angle", angle))
		{
			std::cout << "Warning: Failed to obtain input for angle." << std::endl;
			std::cout << "Using vector: " << angle << std::endl;
		}

		// Check if we have an InstantPulse
		if(this->type != PulseType::InstantPulse)
		{
			// Get pulse time
			if(!this->properties->Get("pulsetime", pulsetime))
			{
				std::cout << "Warning: Failed to obtain input for pulsetime." << std::endl;
				std::cout << "Using vector: " << pulsetime << std::endl;
				
			}

			// Get Pulse field
			if(!this->properties->Get("field", field))
			{
				std::cout << "Warning: Failed to obtain input for the field." << std::endl;
				std::cout << "Using vector: " << field << std::endl;
			}
		}
		else
		{
			std::cout << "Using an Instantpulse; no pulsetime or amplitude is required. " << std::endl;
		}
	}
	
	Pulse::Pulse(const Pulse& _pulse)	: properties(_pulse.properties), type(_pulse.type), group(_pulse.group), rotationaxis(_pulse.rotationaxis), angle(_pulse.angle), pulsetime(_pulse.pulsetime), field(_pulse.field)
	{
	}
	
	Pulse::~Pulse()
	{
	}
	// -----------------------------------------------------
	// Operators
	// -----------------------------------------------------
	const Pulse& Pulse::operator=(const Pulse& _pulse)
	{
		this->properties = std::make_shared<MSDParser::ObjectParser>(*(_pulse.properties));
        this->type = _pulse.type;
        this->group = _pulse.group;
        this->rotationaxis = _pulse.rotationaxis;
        this->angle = _pulse.angle;
        this->pulsetime = _pulse.pulsetime;
        this->field = _pulse.field;
		
		return (*this);
	}
	// ----------------------------------------------------- 
	// Public methods
	// -----------------------------------------------------
	std::string Pulse::Name() const
	{
		return this->properties->Name();
	}
	
	bool Pulse::IsValid() const
	{
		if(this->type == PulseType::Unspecified && !this->group.empty())
			return true;
		else if(this->type == PulseType::InstantPulse && !this->group.empty())
			return true;
        else if(this->type == PulseType::LongPulse && !this->group.empty())
			return true;    
		else if(this->type == PulseType::ShapedPulse && !this->group.empty())
			return true;
		
		return false;
	}
	// -----------------------------------------------------
	// Property methods
	// -----------------------------------------------------
	// Returns the rotation axis
	const arma::vec Pulse::Rotationaxis() const
	{
		return this->rotationaxis;
	}
	
	// Returns the rotation angle 
	const double Pulse::Angle() const
	{
		return this->angle;
	}
	
	// Returns the time of the pulse
	const double Pulse::Pulsetime() const
	{
		return this->pulsetime;
	}

	// Returns the amplitude of the pulse
	const arma::vec Pulse::Field() const
	{
		return this->field;
	}
	
	// Returns pulse type
	PulseType Type(const Pulse& _pulse)
	{
		return _pulse.Type();
	}
	// -----------------------------------------------------
	// Access to custom properties
	// -----------------------------------------------------
	std::shared_ptr<const MSDParser::ObjectParser> Pulse::Properties() const
	{
		return this->properties;
	}
	
	// -----------------------------------------------------
	// Public method to read spins into group,
	// returns false if some spins were not found
	// -----------------------------------------------------
	bool Pulse::ParseSpinGroups(const std::vector<spin_ptr>& _spinlist)
	{
		bool createdSpinLists = false;
		
		// Attempt to get a list of spins from the input file
		std::string str;
		if(!this->properties->Get("spins",str) && !this->properties->Get("group",str))
			return false;
			
		createdSpinLists = this->AddSpinList(str, _spinlist, this->group);
		
		// Check whether we were successful
		if(!createdSpinLists)
		{
			this->group.clear();
			return false;
		}
		
		return true;
	}
	
	// Splits the string into spin names, find the corresponding spins in _spinlist and adds them to _group
	// Private helper method for the ParseSpinGroups method
	bool Pulse::AddSpinList(const std::string& _names, const std::vector<spin_ptr>& _spinlist, std::vector<spin_ptr>& _group, const std::vector<spin_ptr>* _crossCheck)
	{
		// Iterate through the spin name (while-loop) with a ',' delimiter
		std::stringstream ss(_names);
		bool spin_was_found;
		std::string str;
		while(std::getline(ss, str, ','))
		{
			spin_was_found = false;
			
			// Get the spin with the name held by str
			for(auto i = _spinlist.cbegin(); i != _spinlist.cend(); i++)
			{
				if((*i)->Name().compare(str) == 0)
				{
					// Perform crosscheck if requested
					if(_crossCheck != nullptr && std::find((*_crossCheck).cbegin(), (*_crossCheck).cend(), (*i)) != (*_crossCheck).cend())
						return false;
					
					// Make sure we don't have duplicates in the spin group
					if(std::find(_group.cbegin(), _group.cend(), (*i)) != _group.cend())
						return false;
					
					// Add it to the group
					_group.push_back(*i);
					spin_was_found = true;
					break;
				}
			}
			
			if(!spin_was_found)
				return false;
		}
		
		return true;
	}
	// -----------------------------------------------------
	// Public spin list methods
	// -----------------------------------------------------
	
	// -----------------------------------------------------
	// Methods to create ActionTarget objects
	// -----------------------------------------------------
	// Create ActionVectors
	std::vector<RunSection::NamedActionVector> Pulse::CreateActionVectors(const std::string& _system)
	{
		std::vector<RunSection::NamedActionVector> vectors;
		
		if(this->IsValid())
		{

		}
		
		return vectors;
	}
	
	// Create ActionScalars
	std::vector<RunSection::NamedActionScalar> Pulse::CreateActionScalars(const std::string& _system)
	{
		std::vector<RunSection::NamedActionScalar> scalars;
		
		if(this->IsValid())
		{
            // Must be still constructed!
		}
		
		return scalars;
	}
	
	// Method that calls the methods to generate ActionVectors and ActionScalars and inserts them into the given collections
	void Pulse::GetActionTargets(std::vector<RunSection::NamedActionScalar>& _scalars, std::vector<RunSection::NamedActionVector>& _vectors, const std::string& _system)
	{
		// Get ActionTargets from private methods
		auto scalars = this->CreateActionScalars(_system);
		auto vectors = this->CreateActionVectors(_system);
		
		// Insert them
		_scalars.insert(_scalars.end(), scalars.begin(), scalars.end());
		_vectors.insert(_vectors.end(), vectors.begin(), vectors.end());
	}
	// -----------------------------------------------------
}