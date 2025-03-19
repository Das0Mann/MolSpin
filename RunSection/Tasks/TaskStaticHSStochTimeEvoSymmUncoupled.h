/////////////////////////////////////////////////////////////////////////
// TaskStaticHSStochTimeEvoSymmUncoupled (RunSection module) by Gediminas Pazera and Luca Gerhards
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticHSStochTimeEvoSymmUncoupled
#define MOD_RunSection_TaskStaticHSStochTimeEvoSymmUncoupled

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticHSStochTimeEvoSymmUncoupled : public BasicTask
	{
	private:

		double timestep;
		double totaltime;

		SpinAPI::ReactionOperatorType reactionOperators;
		bool productYieldsOnly;			  // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
										  // If false, a quantum yield will be calculated each defined State object
		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticHSStochTimeEvoSymmUncoupled(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticHSStochTimeEvoSymmUncoupled();													// Destructor
	};
}

#endif
