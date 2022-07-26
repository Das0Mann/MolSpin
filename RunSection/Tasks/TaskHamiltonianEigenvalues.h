/////////////////////////////////////////////////////////////////////////
// TaskHamiltonianEigenvalues (RunSection module)
// ------------------
// Prints the eigenvalues, and eigenvectors & Hamiltonian if asked for it.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskHamiltonianEigenvalues
#define MOD_RunSection_TaskHamiltonianEigenvalues

#include "SpinSpace.h"
#include "BasicTask.h"

namespace RunSection
{
	class TaskHamiltonianEigenvalues : public BasicTask
	{
		private:
			// Data members
			bool printEigenvectors;
			bool printHamiltonian;
			bool useSuperspace;
			bool separateRealImag;
			double initialTime;
			double totalTime;
			double timestep;
			bool resonanceFrequencies;
			std::vector<SpinAPI::state_ptr> referenceStates;
			std::vector<std::string> transitionSpins;		// Used for transition matrix element calculations
			
			// Private methods
			void WriteHeader(std::ostream&);	// Write header for the output file
			void GetResonanceFrequencies(const arma::vec&, const arma::cx_mat&, const std::shared_ptr<SpinAPI::SpinSystem>&, const SpinAPI::SpinSpace&);
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskHamiltonianEigenvalues(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskHamiltonianEigenvalues();													// Destructor
	};
}

#endif
