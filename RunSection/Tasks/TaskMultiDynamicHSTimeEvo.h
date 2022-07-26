/////////////////////////////////////////////////////////////////////////
// TaskMultiDynamicHSTimeEvo (RunSection module)
// ------------------
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskMultiDynamicHSTimeEvo
#define MOD_RunSection_TaskMultiDynamicHSTimeEvo

#include "SpinSpace.h"
#include "BasicTask.h"

namespace RunSection
{
	class TaskMultiDynamicHSTimeEvo : public BasicTask
	{
		private:
			double timestep;
			double totaltime;
			unsigned int outputstride;
			
			// Method to obtain creation operators
			void GetCreationOperators(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>>&,
										std::vector<arma::sp_cx_mat>&, const std::vector<arma::cx_mat>&);
			
			// Timestep function
			void AdvanceStep_AsyncLeapfrog(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>>&,
											const std::vector<arma::cx_mat>&, const std::vector<arma::sp_cx_mat>&,
											const std::vector<arma::cx_mat>&, const std::vector<arma::sp_cx_mat>&,
											const std::vector<arma::sp_cx_mat>&, std::vector<arma::cx_mat>&);
			
			// Method to obtain the results from the current state
			void OutputResults(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>>&, const std::vector<arma::cx_mat>&, const unsigned int);
			
			// Method to update time-dependent interactions or reactions
			void UpdateTimeDependences(const std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>>&,
										std::vector<arma::sp_cx_mat>&, std::vector<arma::sp_cx_mat>&, unsigned int);
			
			// Write header for the output file
			void WriteHeader(std::ostream&);
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskMultiDynamicHSTimeEvo(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskMultiDynamicHSTimeEvo();													// Destructor
	};
}

#endif
