/////////////////////////////////////////////////////////////////////////
// TaskStaticSSRelaxation (RunSection module)
// ------------------
// 
// Simple quantum yield calculation in Liouville space, derived from the
// properties of the Laplace transformation.
// 
// (c) 2018 by Claus N.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSRelaxation
#define MOD_RunSection_TaskStaticSSRelaxation

#include "BasicTask.h"

namespace RunSection
{
	class TaskStaticSSRelaxation : public BasicTask
	{
		private:
			double relaxationRate;
			std::vector<std::string> relaxingSpins;
			bool productYieldsOnly;				// If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
												// If false, a quantum yield will be calculated each defined State object
			
			void WriteHeader(std::ostream&);	// Write header for the output file
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskStaticSSRelaxation(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskStaticSSRelaxation();													// Destructor
	};
}

#endif
