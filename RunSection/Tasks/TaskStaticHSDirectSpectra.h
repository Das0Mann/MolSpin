/////////////////////////////////////////////////////////////////////////
// TaskStaticHSDirectSpectra (RunSection module) by Luca Gerhards
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticHSDirectSpectra
#define MOD_RunSection_TaskStaticHSDirectSpectra

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticHSDirectSpectra : public BasicTask
	{
	private:

		double timestep;
		double totaltime;

		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticHSDirectSpectra(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticHSDirectSpectra();												   // Destructor
	};
}

#endif
