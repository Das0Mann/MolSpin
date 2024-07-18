/////////////////////////////////////////////////////////////////////////
// TaskStaticHSSymmetricDecay (RunSection module)
// ------------------
// Hilbert space method. Assummes symmetric (spin independent) recombination.
// Note: Also assumes that we have a radical pair with non-interacting spins!
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticRPOnlyHSSymDecRedfield
#define MOD_RunSection_TaskStaticRPOnlyHSSymDecRedfield

#include "SpinSpace.h"
#include "BasicTask.h"

namespace RunSection
{
	class TaskStaticRPOnlyHSSymDecRedfield : public BasicTask
	{
	private:
		bool modeQuantumYield;
		bool productYieldsOnly; // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
		double timestep;
		double totaltime;
		// If false, a quantum yield will be calculated each defined State object

		void WriteHeader(std::ostream &); // Write header for the output file
		void GetQuantumYields(const SpinAPI::SpinSpace &_space, const arma::cx_mat &_rho, const unsigned int _step, const std::vector<SpinAPI::state_ptr> &_states, const std::vector<SpinAPI::transition_ptr> &_transitions, std::map<std::string, arma::vec> &_yields);
		void OutputQuantumYields(const SpinAPI::SpinSpace &_space, std::map<std::string, arma::vec> &_yields, const unsigned int _steps, const std::vector<SpinAPI::state_ptr> &_states, const std::vector<SpinAPI::transition_ptr> &_transitions);

		bool Redfieldtensor(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, arma::cx_mat &_redfieldtensor);
		bool ConstructSpecDensGeneral(const int &_spectral_function, const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_domega, arma::cx_mat &_specdens);
		bool ConstructSpecDensSpecific(const int &_spectral_function, const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_domega, arma::cx_mat &_specdens);
		bool Slippage(arma::cx_mat **_ptr_Tensors, const int &_num_op, const arma::cx_mat &_eig_val_mat, const arma::cx_mat &_domega, const arma::cx_mat &_rho0, const std::complex<double> &_tau_c, arma::cx_mat &_rho0_new);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticRPOnlyHSSymDecRedfield(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticRPOnlyHSSymDecRedfield();												   // Destructor
	};
}

#endif
