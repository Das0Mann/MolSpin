/////////////////////////////////////////////////////////////////////////
// TaskMultiStaticSSNakajimaZwanzigTimeEvo (RunSection module)
// ------------------
//
// Simple time-evolution calculation in Liouville space.
//
// -- Multi-system version: Allows transitions between SpinSystems --
//
// Molecular Spin Dynamics Software - developed by Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiStaticSSNakajimaZwanzigTimeEvo
#define MOD_RunSection_TaskMultiStaticSSNakajimaZwanzigTimeEvo

#include "BasicTask.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskMultiStaticSSNakajimaZwanzigTimeEvo : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

		// Private method that gathers and outputs the results from a given time-integrated density operator
		bool NakajimaZwanzigtensorTimeEvo(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, arma::cx_mat &_NakajimaZwanzigtensor); // Contruction of NakajimaZwanzigtensor with operator basis
		bool ConstructSpecDensGeneralTimeEvo(const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_omega, arma::cx_mat &_specdens);
		bool ConstructSpecDensSpecificTimeEvo(const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_omega, arma::cx_mat &_specdens);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskMultiStaticSSNakajimaZwanzigTimeEvo(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskMultiStaticSSNakajimaZwanzigTimeEvo();												   // Destructor
	};
}

#endif
