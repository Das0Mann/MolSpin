/////////////////////////////////////////////////////////////////////////
// TaskStaticSSSpectraNakajimaZwanzig (RunSection module)  developed by Irina Anisimova and Luca Gerhards.
// ------------------
//
// Simple quantum yield calculation in Liouville space, derived from the
// properties of the Laplace transformation.
//
// Molecular Spin Dynamics Software - developed by Irina Anisimova and Luca Gerhards.
// (c) 2022 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSSpectraNakajimaZwanzig
#define MOD_RunSection_TaskStaticSSSpectraNakajimaZwanzig

#include "BasicTask.h"
#include "SpinAPIDefines.h"
#include "SpinSpace.h"

namespace RunSection
{
	//Because the declaration of i is long define it as a new variable that is easier to use
	using SystemIterator = std::vector<SpinAPI::system_ptr>::const_iterator;

	class TaskStaticSSSpectraNakajimaZwanzig : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

		bool NakajimaZwanzigtensorSpectra(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, arma::cx_mat &_NakajimaZwanzigtensor); // Contruction of NakajimaZwanzigtensor with operator basis
		bool ConstructSpecDensGeneralSpectra(const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_omega, arma::cx_mat &_specdens);
		bool ConstructSpecDensSpecificSpectra(const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_omega, arma::cx_mat &_specdens);

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSSpectraNakajimaZwanzig(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSSpectraNakajimaZwanzig();													// Destructor
	};

}

#endif
