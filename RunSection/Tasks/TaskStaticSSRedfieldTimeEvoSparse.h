//////////////////////////////////////////////////////////////////////////////
// MolSpin - Redfield Theory Task - developed by Luca Gerhards
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSRedfieldTimeEvoSparse
#define MOD_RunSection_TaskStaticSSRedfieldTimeEvoSparse

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticSSRedfieldTimeEvoSparse : public BasicTask
	{
		private:
			double timestep;
			double totaltime;
			SpinAPI::ReactionOperatorType reactionOperators;
			
			void WriteHeader(std::ostream&);	// Write header for the output file
			bool RedfieldtensorTimeEvoSparse(const arma::sp_cx_mat& _op1, const arma::sp_cx_mat& _op2, const arma::sp_cx_mat& _specdens,arma::sp_cx_mat& _redfieldtensor); //Contruction of Redfieldtensor with operator basis
			bool ConstructSpecDensGeneralTimeEvoSparse(const int& _spectral_function, const std::vector<double>& _ampl_list, const std::vector<double>& _tau_c_list,const arma::sp_cx_mat& _domega, arma::sp_cx_mat& _specdens);
			bool ConstructSpecDensSpecificTimeEvoSparse(const int& _spectral_function, const std::complex<double>& _ampl, const std::complex<double>& _tau_c,const arma::sp_cx_mat& _domega, arma::sp_cx_mat& _specdens);

			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskStaticSSRedfieldTimeEvoSparse(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskStaticSSRedfieldTimeEvoSparse();													// Destructor
	};
}

#endif
