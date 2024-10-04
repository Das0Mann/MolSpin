/////////////////////////////////////////////////////////////////////////
// ActionAddVector implementation (RunSection module)
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include "ActionAddVector.h"
#include "ObjectParser.h"

namespace RunSection
{
	// -----------------------------------------------------
	// ActionAddVector Constructors and Destructor
	// -----------------------------------------------------
	ActionAddVector::ActionAddVector(const MSDParser::ObjectParser &_parser, const std::map<std::string, ActionScalar> &_scalars, const std::map<std::string, ActionVector> &_vectors) : Action(_parser, _scalars, _vectors), actionVector(nullptr), direction()
	{
	}

	ActionAddVector::~ActionAddVector()
	{
	}
	// -----------------------------------------------------
	// ActionAddVector derived methods
	// -----------------------------------------------------
	// Method to perform a step
	bool ActionAddVector::DoStep()
	{
		// Make sure we have an ActionVector to act on
		if (actionVector == nullptr || !this->IsValid())
		{
			return false;
		}

		// Retrieve the vector we want to change
		arma::vec vec = actionVector->Get();

		// Add a vector to the ActionVector
		vec += this->Value() * this->direction;

		// Set the new vector
		return this->actionVector->Set(vec);
	}

	// Method to prepare the action and check whether it is valid
	bool ActionAddVector::DoValidate()
	{
		// Make sure that a direction is specified
		arma::vec tmp;
		if (!this->Properties()->Get("direction", tmp))
		{
			std::cout << "ERROR: No direction specified for the AddVector action \"" << this->Name() << "\"!" << std::endl;
			return false;
		}

		// Attemp to set the direction
		if (!this->SetDirection(tmp))
		{
			std::cout << "ERROR: Invalid direction specified for the AddVector action \"" << this->Name() << "\"!" << std::endl;
			return false;
		}

		// Get the ActionVector name
		std::string str;
		if (!this->Properties()->Get("actionvector", str) && !this->Properties()->Get("vector", str))
		{
			std::cout << "ERROR: No ActionVector specified for the AddVector action \"" << this->Name() << "\"!" << std::endl;
			return false;
		}

		// Attemp to set the ActionVector
		if (!this->Vector(str, &(this->actionVector)))
		{
			std::cout << "ERROR: Could not find ActionVector \"" << str << "\" specified for the AddVector action \"" << this->Name() << "\"!" << std::endl;
			return false;
		}

		// Readonly ActionTargets cannot be acted on
		if (this->actionVector->IsReadonly())
		{
			std::cout << "ERROR: Readonly ActionVector \"" << str << "\" specified for the AddVector action \"" << this->Name() << "\"! Cannot act on this vector!" << std::endl;
			return false;
		}

		return true;
	}
	// -----------------------------------------------------
	// ActionAddVector private methods
	// -----------------------------------------------------
	// Sets the direction unit vector
	bool ActionAddVector::SetDirection(const arma::vec &_vec)
	{
		// Check whether we have a valid direction vector
		if (_vec.n_elem == 3 && _vec.is_finite())
		{
			// Set the direction, and normalize it (the length is set by "value" on the action)
			this->direction = normalise(_vec);

			return true;
		}

		return false;
	}

	bool ActionAddVector::Reset()
	{
		this->actionVector->Reset();
	}
	// -----------------------------------------------------
}
