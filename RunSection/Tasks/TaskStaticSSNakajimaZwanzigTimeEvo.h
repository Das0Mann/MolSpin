//////////////////////////////////////////////////////////////////////////////
// MolSpin - Nakajima Zwanzig Theory Time Evolution Task - developed by Luca Gerhards
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSNakajimaZwanzigTimeEvo
#define MOD_RunSection_TaskStaticSSNakajimaZwanzigTimeEvo

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticSSNakajimaZwanzigTimeEvo : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &);																																						  // Write header for the output file
		bool NakajimaZwanzigtensorTimeEvo(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, arma::cx_mat &_NakajimaZwanzigtensor); // Contruction of NakajimaZwanzigtensor with operator basis
		bool ConstructSpecDensGeneralTimeEvo(const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_omega, arma::cx_mat &_specdens);
		bool ConstructSpecDensSpecificTimeEvo(const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_omega, arma::cx_mat &_specdens);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSNakajimaZwanzigTimeEvo(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSNakajimaZwanzigTimeEvo();													 // Destructor
	};
}

#endif