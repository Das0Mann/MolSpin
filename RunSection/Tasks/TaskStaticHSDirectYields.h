/////////////////////////////////////////////////////////////////////////
// TaskStaticHSDirectYields (RunSection module) by Gediminas Pazera
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticHSDirectYields
#define MOD_RunSection_TaskStaticHSDirectYields

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticHSDirectYields : public BasicTask
	{
		private:
			SpinAPI::ReactionOperatorType reactionOperators;
			bool productYieldsOnly;				// If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
												// If false, a quantum yield will be calculated each defined State object
			void WriteHeader(std::ostream&);	// Write header for the output file
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskStaticHSDirectYields(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskStaticHSDirectYields();													// Destructor
	};
}

#endif
