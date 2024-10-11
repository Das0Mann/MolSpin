/////////////////////////////////////////////////////////////////////////
// TaskStaticSSSpectra (RunSection module)  developed by Irina Anisimova and Luca Gerhards.
// ------------------
//
// Simple quantum yield calculation in Liouville space, derived from the
// properties of the Laplace transformation.
//
// Molecular Spin Dynamics Software - developed by Irina Anisimova and Luca Gerhards.
// (c) 2022 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_TaskStaticSSSpectra
#define MOD_RunSection_TaskStaticSSSpectra

#include "BasicTask.h"
#include "SpinAPIDefines.h"
#include "SpinSpace.h"

namespace RunSection
{
	//Because the declaration of i is long define it as a new variable that is easier to use
	using SystemIterator = std::vector<SpinAPI::system_ptr>::const_iterator;

	class TaskStaticSSSpectra : public BasicTask
	{
	private:
		double timestep;
		double totaltime;
		SpinAPI::ReactionOperatorType reactionOperators;

		void WriteHeader(std::ostream &); // Write header for the output file

	protected:
		bool RunLocal() override;
		bool Validate() override;

	public:
		// Constructors / Destructors
		TaskStaticSSSpectra(const MSDParser::ObjectParser &, const RunSection &); // Normal constructor
		~TaskStaticSSSpectra();													// Destructor
	};

}

#endif
