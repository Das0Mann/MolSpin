/////////////////////////////////////////////////////////////////////////
// TaskMultiStaticSSRedfieldTimeEvo (RunSection module)
// ------------------
//
// Simple time-evolution calculation in Liouville space.
//
// -- Multi-system version: Allows transitions between SpinSystems --
//
// Molecular Spin Dynamics Software - developed by Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiStaticSSRedfieldTimeEvo
#define MOD_RunSection_TaskMultiStaticSSRedfieldTimeEvo

#include "BasicTask.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskMultiStaticSSRedfieldTimeEvo : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

		// Private method that gathers and outputs the results from a given time-integrated density operator
		void GatherResults(const arma::cx_mat &, const SpinAPI::SpinSystem &, const SpinAPI::SpinSpace &);
		bool Redfieldtensor(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, arma::cx_mat &_redfieldtensor);															   // Contruction of Redfieldtensor with operator basis
		bool ConstructSpecDensGeneral(const int &_spectral_function, const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_domega, arma::cx_mat &_specdens); // Construction of Spectral Density
		bool ConstructSpecDensSpecific(const int &_spectral_function, const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_domega, arma::cx_mat &_specdens);
		bool Slippage(arma::cx_mat **_ptr_Tensors, const int &_num_op, const arma::cx_mat &_eig_val_mat, const arma::cx_mat &_domega, const arma::cx_mat &_rho0, const std::complex<double> &_tau_c, arma::cx_mat &_rho0_new);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskMultiStaticSSRedfieldTimeEvo(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskMultiStaticSSRedfieldTimeEvo();												   // Destructor
	};
}

#endif
