/////////////////////////////////////////////////////////////////////////
// TaskMultiStaticSSTimeEvoSpectra (RunSection module)
// ------------------
//
// Simple time-evolution calculation in Liouville space.
//
// -- Multi-system version: Allows transitions between SpinSystems --
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiStaticSSTimeEvoSpectra
#define MOD_RunSection_TaskMultiStaticSSTimeEvoSpectra

#include "BasicTask.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskMultiStaticSSTimeEvoSpectra : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		bool productYieldsOnly;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

		// Private method that gathers and outputs the results from a given time-integrated density operator
		bool GatherResults(const arma::cx_mat &, const SpinAPI::SpinSystem &, const SpinAPI::SpinSpace &);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskMultiStaticSSTimeEvoSpectra(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskMultiStaticSSTimeEvoSpectra();												   // Destructor
	};
}

#endif
