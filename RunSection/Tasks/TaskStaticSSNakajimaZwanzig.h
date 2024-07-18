//////////////////////////////////////////////////////////////////////////////
// MolSpin - Nakajima Zwanzig Theory Task - developed by Luca Gerhards
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSNakajimaZwanzig
#define MOD_RunSection_TaskStaticSSNakajimaZwanzig

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticSSNakajimaZwanzig : public BasicTask
	{
	private:
		SpinAPI::ReactionOperatorType reactionOperators;
		bool productYieldsOnly;																																							   // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
																																														   // If false, a quantum yield will be calculated each defined State object
		void WriteHeader(std::ostream &);																																				   // Write header for the output file
		bool NakajimaZwanzigtensor(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, const arma::cx_mat _eigenvec, arma::cx_mat &_NakajimaZwanzigtensor); // Contruction of NakajimaZwanzigtensor with operator basis
		bool ConstructSpecDensGeneral(const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_omega, arma::cx_mat &_specdens);
		bool ConstructSpecDensSpecific(const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_omega, arma::cx_mat &_specdens);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSNakajimaZwanzig(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSNakajimaZwanzig();													  // Destructor
	};
}

#endif