/////////////////////////////////////////////////////////////////////////
// TaskMultiRadicalPairSSTimeEvo (RunSection module)
// ------------------
//
// Simple time-evolution calculation in Liouville space.
//
// -- Multi-system version: Allows transitions between SpinSystems --
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiRadicalPairSSTimeEvo
#define MOD_RunSection_TaskMultiRadicalPairSSTimeEvo

#include "BasicTask.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskMultiRadicalPairSSTimeEvo : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

		// Private method that gathers and outputs the results from a given time-integrated density operator
		void GatherResults(const arma::cx_mat &, const SpinAPI::SpinSystem &, const SpinAPI::SpinSpace &);
		bool GenerateHamiltonian(const std::vector<SpinAPI::interaction_ptr> interactions, arma::sp_cx_mat& H, int dimension, std::shared_ptr<SpinAPI::SpinSpace> SpinSystem);
		bool GenerateReactionOperator(const std::vector<SpinAPI::transition_ptr> transitions, arma::sp_cx_mat& K, int dimension, std::shared_ptr<SpinAPI::SpinSpace> SpinSystem);
	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskMultiRadicalPairSSTimeEvo(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskMultiRadicalPairSSTimeEvo();												   // Destructor
	};
}

#endif
