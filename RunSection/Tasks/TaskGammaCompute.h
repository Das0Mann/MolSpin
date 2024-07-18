/////////////////////////////////////////////////////////////////////////
// TaskGammaCompute (RunSection module)
// ------------------
// An implementation of the Gamma-COMPUTE algorithm for calculating RF
// magnetic field effects with phase averaging.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskGammaCompute
#define MOD_RunSection_TaskGammaCompute

#include "SpinSpace.h"
#include "BasicTask.h"

namespace RunSection
{
	class TaskGammaCompute : public BasicTask
	{
	private:
		unsigned int steps;
		double totaltime;

		// Private methods
		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskGammaCompute(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskGammaCompute();												   // Destructor
	};
}

#endif
