/////////////////////////////////////////////////////////////////////////
// ActionAddScalar (RunSection module)
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ActionAddScalar
#define MOD_RunSection_ActionAddScalar

#include "Action.h"

namespace RunSection
{
	class ActionAddScalar : public Action
	{
		private:
			// Data members
			ActionScalar* actionScalar;
			
		protected:
			// Overwritten protected methods
			bool DoStep() override;
			bool DoValidate() override;
			
		public:
			// Constructors / Destructors
			ActionAddScalar(const MSDParser::ObjectParser&, const std::map<std::string, ActionScalar>&, const std::map<std::string, ActionVector>&);		// Normal constructor
			~ActionAddScalar();
	};
}

#endif
