/////////////////////////////////////////////////////////////////////////
// TaskPeriodicHSTimeEvo (RunSection module)
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskPeriodicHSTimeEvo
#define MOD_RunSection_TaskPeriodicHSTimeEvo

#include "SpinSpace.h"
#include "BasicTask.h"

namespace RunSection
{
	class TaskPeriodicHSTimeEvo : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		unsigned int stepsPerPeriod;
		unsigned int outputstride;
		bool modeQuantumYield;
		bool productYieldsOnly; // If true, a quantum yield will be calculated from each Transition object and multiplied by the rate constant
								// If false, a quantum yield will be calculated each defined State object

		bool timedependentInteractions;

		// Timestep function
		void AdvanceStep_AsyncLeapfrog(const arma::cx_mat &, const arma::cx_mat &, arma::cx_mat &);
		void OutputTimeEvolution(const SpinAPI::SpinSpace &, const arma::cx_mat &, const unsigned int, const std::vector<SpinAPI::state_ptr> &);
		void GetQuantumYields(const SpinAPI::SpinSpace &, const arma::cx_mat &, const unsigned int, const std::vector<SpinAPI::state_ptr> &, const std::vector<SpinAPI::transition_ptr> &, std::map<std::string, arma::vec> &);
		void OutputQuantumYields(const SpinAPI::SpinSpace &, std::map<std::string, arma::vec> &, const unsigned int, const std::vector<SpinAPI::state_ptr> &, const std::vector<SpinAPI::transition_ptr> &);
		double GetPeriod(const SpinAPI::system_ptr &);
		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskPeriodicHSTimeEvo(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskPeriodicHSTimeEvo();													// Destructor
	};
}

#endif
