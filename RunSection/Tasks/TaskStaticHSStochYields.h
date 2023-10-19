/////////////////////////////////////////////////////////////////////////
// TaskStaticHSStochYields (RunSection module) by Gediminas Pazera
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticHSStochYields
#define MOD_RunSection_TaskStaticHSStochYields

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticHSStochYields : public BasicTask
	{
		private:
			SpinAPI::ReactionOperatorType reactionOperators;
			bool productYieldsOnly;				// If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
												// If false, a quantum yield will be calculated each defined State object
		        arma::cx_colvec KrylovExpmGeneral(const arma::sp_cx_mat &H, const arma::cx_colvec &b, const arma::cx_double dt, int KryDim, int HilbSize);
			arma::cx_colvec KrylovExpmSymm(const arma::sp_cx_mat &H, const arma::cx_colvec &b, const arma::cx_double dt, int KryDim, int HilbSize); 	
			void ArnoldiProcess(const arma::sp_cx_mat &H, const arma::cx_colvec &b, arma::cx_mat &KryBasis, arma::cx_mat &Hessen, int KryDim, double& h_mplusone_m); // Produce Krylov Basis by Arnoldi Process
			void LanczosProcess(const arma::sp_cx_mat &H, const arma::cx_colvec &b, arma::cx_mat &KryBasis, arma::cx_mat &Hessen, int KryDim, double& h_mplusone_m); // Produce Krylov Basis by Lanczos Process
			void WriteHeader(std::ostream&);	// Write header for the output file
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskStaticHSStochYields(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskStaticHSStochYields();													// Destructor
	};
}

#endif
