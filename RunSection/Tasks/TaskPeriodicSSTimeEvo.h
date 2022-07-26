/////////////////////////////////////////////////////////////////////////
// TaskPeriodicSSTimeEvo (RunSection module)
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskPeriodicSSTimeEvo
#define MOD_RunSection_TaskPeriodicSSTimeEvo

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskPeriodicSSTimeEvo : public BasicTask
	{
		private:
			double timestep;
			double totaltime;
			SpinAPI::ReactionOperatorType reactionOperators;
			unsigned int stepsPerPeriod;
			
			void WriteHeader(std::ostream&);	// Write header for the output file
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskPeriodicSSTimeEvo(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskPeriodicSSTimeEvo();													// Destructor
	};
}

#endif
