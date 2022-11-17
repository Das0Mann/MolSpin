//////////////////////////////////////////////////////////////////////////////
// MolSpin - Redfield Theory Task - developed by Luca Gerhards
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSSpectra
#define MOD_RunSection_TaskStaticSSSpectra

#include "BasicTask.h"
#include "SpinAPIDefines.h"

namespace RunSection
{
	class TaskStaticSSSpectra : public BasicTask
	{
		private:
			double timestep;
			double totaltime;
			SpinAPI::ReactionOperatorType reactionOperators;
			
			void WriteHeader(std::ostream&);	// Write header for the output file
			
		protected:
			bool RunLocal() override;
			bool Validate() override;
			
		public:
			// Constructors / Destructors
			TaskStaticSSSpectra(const MSDParser::ObjectParser&, const RunSection&);	// Normal constructor
			~TaskStaticSSSpectra();													// Destructor
	};
}

#endif
