/////////////////////////////////////////////////////////////////////////
// ActionMultiplyScalar (RunSection module)
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ActionMultiplyScalar
#define MOD_RunSection_ActionMultiplyScalar

#include "Action.h"

namespace RunSection
{
	class ActionMultiplyScalar : public Action
	{
	private:
		// Data members
		ActionScalar *actionScalar;

	protected:
		// Overwritten protected methods
		bool DoStep() override;
		bool DoValidate() override;
		bool Reset() override;

	public:
		// Constructors / Destructors
		ActionMultiplyScalar(const MSDParser::ObjectParser &, const std::map<std::string, ActionScalar> &, const std::map<std::string, ActionVector> &); // Normal constructor
		~ActionMultiplyScalar();
	};
}

#endif
