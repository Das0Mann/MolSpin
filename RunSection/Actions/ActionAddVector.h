/////////////////////////////////////////////////////////////////////////
// ActionAddVector (RunSection module)
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ActionAddVector
#define MOD_RunSection_ActionAddVector

#include "Action.h"

namespace RunSection
{
	class ActionAddVector : public Action
	{
	private:
		// Data members
		ActionVector *actionVector;
		arma::vec direction;

		// Private methods
		bool SetDirection(const arma::vec &);

	protected:
		// Overwritten protected methods
		bool DoStep() override;
		bool DoValidate() override;
		bool Reset() override;

	public:
		// Constructors / Destructors
		ActionAddVector(const MSDParser::ObjectParser &, const std::map<std::string, ActionScalar> &, const std::map<std::string, ActionVector> &); // Normal constructor
		~ActionAddVector();
	};
}

#endif
