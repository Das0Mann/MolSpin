//////////////////////////////////////////////////////////////////////////////
// MolSpin - Redfield Theory Task - developed by Luca Gerhards
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSRedfieldSparse
#define MOD_RunSection_TaskStaticSSRedfieldSparse

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticSSRedfieldSparse : public BasicTask
	{
	private:
		SpinAPI::ReactionOperatorType reactionOperators;
		bool productYieldsOnly; // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
								// If false, a quantum yield will be calculated each defined State object

		void WriteHeader(std::ostream &);																														 // Write header for the output file
		bool RedfieldtensorSparse(const arma::sp_cx_mat &_op1, const arma::sp_cx_mat &_op2, const arma::sp_cx_mat &_specdens, arma::sp_cx_mat &_redfieldtensor); // Contruction of Redfieldtensor with operator basis
		bool ConstructSpecDensGeneralSparse(const int &_spectral_function, const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::sp_cx_mat &_domega, arma::sp_cx_mat &_specdens);
		bool ConstructSpecDensSpecificSparse(const int &_spectral_function, const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::sp_cx_mat &_domega, arma::sp_cx_mat &_specdens);
		bool SlippageSparse(arma::sp_cx_mat **_ptr_Tensors, const int &_num_op, const arma::sp_cx_mat &_eig_val_mat, const arma::sp_cx_mat &_domega, const arma::sp_cx_mat &_rho0, const std::complex<double> &_tau_c, arma::sp_cx_mat &_rho0_new);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSRedfieldSparse(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSRedfieldSparse();													 // Destructor
	};
}

#endif
