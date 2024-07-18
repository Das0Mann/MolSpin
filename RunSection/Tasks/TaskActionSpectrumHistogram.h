/////////////////////////////////////////////////////////////////////////
// TaskActionSpectrumHistogram (RunSection module)
// ------------------
// Computes the action-spectrum histogram of a spin system, see Hiscock et al. (2017) and Leberecht et al. (2022).
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// Task implemented by Siu Ying Wong and Luca Gerhards
// (c) 2022 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskActionSpectrumHistogram
#define MOD_RunSection_TaskActionSpectrumHistogram

#include "SpinSpace.h"
#include "BasicTask.h"

namespace RunSection
{
	class TaskActionSpectrumHistogram : public BasicTask
	{
	private:
		// Data members

		// Private methods
		arma::vec GetHistogramBinCenters();
		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskActionSpectrumHistogram(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskActionSpectrumHistogram();													  // Destructor
	};
}

#endif
