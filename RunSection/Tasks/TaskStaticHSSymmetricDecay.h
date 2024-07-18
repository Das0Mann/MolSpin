/////////////////////////////////////////////////////////////////////////
// TaskStaticHSSymmetricDecay (RunSection module)
// ------------------
// Hilbert space method. Assummes symmetric (spin independent) recombination.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticHSSymmetricDecay
#define MOD_RunSection_TaskStaticHSSymmetricDecay

#include "BasicTask.h"

namespace RunSection
{
	class TaskStaticHSSymmetricDecay : public BasicTask
	{
	private:
		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticHSSymmetricDecay(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticHSSymmetricDecay();													 // Destructor
	};
}

#endif
