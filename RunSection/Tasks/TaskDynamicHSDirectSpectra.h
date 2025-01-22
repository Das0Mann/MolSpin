/////////////////////////////////////////////////////////////////////////
// TaskDynamicHSDirectSpectra (RunSection module) by Luca Gerhards
// ------------------
//
// Very efficient quantum yield calculations in Hilbert Space, using Monte Carlo sampling.
// For dynamic time-dependent Hamiltonians.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskDynamicHSDirectSpectra
#define MOD_RunSection_TaskDynamicHSDirectSpectra

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskDynamicHSDirectSpectra : public BasicTask
	{
	private:
		// Things for time dependent Hamiltonian

		double timestep;
		double totaltime;
		bool timedependentInteractions;
		bool timedependentTransitions;

		//

		SpinAPI::ReactionOperatorType reactionOperators;
		bool productYieldsOnly;			  // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
										  // If false, a quantum yield will be calculated each defined State object
		void WriteHeader(std::ostream &); // Write header for the output file
	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskDynamicHSDirectSpectra(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskDynamicHSDirectSpectra();													// Destructor
		bool is_identity_matrix(arma::sp_cx_mat &matrix);								// check if the matrix is an identity
	};
}

#endif
