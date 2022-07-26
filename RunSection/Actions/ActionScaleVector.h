/////////////////////////////////////////////////////////////////////////
// ActionScaleVector (RunSection module)
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ActionScaleVector
#define MOD_RunSection_ActionScaleVector

#include "Action.h"

namespace RunSection
{
	class ActionScaleVector : public Action
	{
		private:
			// Data members
			ActionVector* actionVector;
			
		protected:
			// Overwritten protected methods
			bool DoStep() override;
			bool DoValidate() override;
			
		public:
			// Constructors / Destructors
			ActionScaleVector(const MSDParser::ObjectParser&, const std::map<std::string, ActionScalar>&, const std::map<std::string, ActionVector>&);		// Normal constructor
			~ActionScaleVector();
	};
}

#endif
