/////////////////////////////////////////////////////////////////////////
// ActionRotateVector (RunSection module)
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_ActionRotateVector
#define MOD_RunSection_ActionRotateVector

#include "Action.h"

namespace RunSection
{
	class ActionRotateVector : public Action
	{
		private:
			// Data members
			ActionVector* actionVector;
			arma::vec axis;
			arma::vec refAxis1;
			arma::vec refAxis2;
			
			// Private methods
			bool SetAxis(const arma::vec&);
			
		protected:
			// Overwritten protected methods
			bool DoStep() override;
			bool DoValidate() override;
			
		public:
			// Constructors / Destructors
			ActionRotateVector(const MSDParser::ObjectParser&, const std::map<std::string, ActionScalar>&, const std::map<std::string, ActionVector>&);		// Normal constructor
			~ActionRotateVector();
	};
}

#endif
