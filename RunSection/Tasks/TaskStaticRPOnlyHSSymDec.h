/////////////////////////////////////////////////////////////////////////
// TaskStaticHSSymmetricDecay (RunSection module)
// ------------------
// Hilbert space method. Assummes symmetric (spin independent) recombination.
// Note: Also assumes that we have a radical pair with non-interacting spins!
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticRPOnlyHSSymDec
#define MOD_RunSection_TaskStaticRPOnlyHSSymDec

#include "BasicTask.h"

namespace RunSection
{
	class TaskStaticRPOnlyHSSymDec : public BasicTask
	{
		private:
			void WriteHeader(std::ostream&);	// Write header for the output file
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskStaticRPOnlyHSSymDec(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskStaticRPOnlyHSSymDec();													// Destructor
	};
}

#endif
