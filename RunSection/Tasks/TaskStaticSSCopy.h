/////////////////////////////////////////////////////////////////////////
// TaskStaticSS (RunSection module)
// ------------------
//
// Simple quantum yield calculation in Liouville space, derived from the
// properties of the Laplace transformation.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSCopy
#define MOD_RunSection_TaskStaticSSCopy

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticSSCopy : public BasicTask
	{
	private:
		SpinAPI::ReactionOperatorType reactionOperators;
		bool productYieldsOnly; // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
								// If false, a quantum yield will be calculated each defined State object

		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSCopy(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSCopy();												   // Destructor
	};
}

#endif
