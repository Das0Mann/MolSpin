/////////////////////////////////////////////////////////////////////////
// MSDParser class (MSDParser Module)
// ------------------
// Class used to parse MSD input files.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_MSDParser_MSDParser
#define MOD_MSDParser_MSDParser

#include <vector>
#include <memory>
#include "MSDParserDefines.h"
#include "SpinAPIfwd.h"
#include "ObjectParser.h"
#include "RunSectionfwd.h"

namespace MSDParser
{
	class MSDParser
	{
		private:
			// Implementation
			bool isLoaded = false;
			std::string filename;
			std::vector<std::shared_ptr<SpinAPI::SpinSystem>> systems;
			std::vector<ObjectParser> runTasks;
			std::vector<ObjectParser> runActions;
			std::vector<ObjectParser> settingsObjects;
			std::vector<SpinAPI::output_ptr> runOutputs;
		
		public:
			// Constructors / Destructors
			explicit MSDParser(std::string);		// Normal constructor
			MSDParser(const MSDParser&) = delete;	// Copy-constructor
			~MSDParser();							// Destructor
			
			// Operators
			MSDParser& operator=(const MSDParser&) = delete;	// Copy-assignment
			
			// Public methods
			bool Load();
			std::string Filename() const {return this->filename;};
			void FillRunSection(RunSection::RunSection&);		// Loads data into a RunSection object
	};
	
	// Non-member non-friend functions
	bool FileExists(const MSDParser&);	// Check whether the file passed to the MSDParser exists
	bool FileLoaded(const MSDParser&);	// Check whether the file passed to the MSDParser has been loaded
}

#endif
