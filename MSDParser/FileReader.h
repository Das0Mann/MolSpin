/////////////////////////////////////////////////////////////////////////
// FileReader class (MSDParser Module)
// ------------------
// Reads files and passes the read objects on to the MSDParser.
// 
// Any errors in reading the files will be handled by the FileReader,
// and if the file could not be found/opened, an object with
// ObjectType::EOF will be returned by ReadObject() and an error message
// will be written by FileReader.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_MSDParser_FileReader
#define MOD_MSDParser_FileReader

#include <map>
#include <vector>
#include <memory>
#include <fstream>
#include "MSDParserDefines.h"

namespace MSDParser
{
	class FileReader
	{
		private:
			// Static list of read files to avoid circular includes
			static std::vector<std::string> filelist;
			static std::map<std::string, std::string> define_directives;
			
			// Static methods
			static std::string PathFromFileName(const std::string&);
			static bool PathIsAbsolute(const std::string&);
			
			// Implementation details
			ObjectGroup objectGroupType;
			std::string objectGroupName;
			std::shared_ptr<FileReader> subfile;
			std::ifstream filehandle;		// Initialized in constructor
			std::string filename;			// Used when reporting errors
			int linenumber;					// Keeps track of linenumber in the file held by the FileReader (used to report parsing errors)
			bool inObjectGroup;				// Whether or not the FileReader file pointer is inside an ObjectGroup
			bool inExcludedDirective;		// Whether or not the FileReader file pointer is inside an "ifdef", "ifndef" or "else" directive that should be skipped (e.g. "ifdef" not satisfied)
			
			// Private methods to validate read characters
			bool isAlphanumericCharacter(const char&) const;
			bool isAllowedCharacter(const char&) const;
			bool isAllowedCharacterInString(const char&) const;
			char getLowerCase(const char&) const;	// Converts uppercase letters to lowercase, works similar to isAlphanumericCharacter
			
			// Private method to handle '#'-directives read from the input file
			void handleDirective(const std::string&,const std::string&,const std::string&);
			
			// Private method to set the current object group
			void setObjectGroup(const std::string&,const std::string&);
			
			// Private method to create an object, once enough data is read
			MSDFileObject createObject(const std::string&,const std::string&,std::string);
		
		public:
			// Constructors / Destructors
			explicit FileReader(const std::string&);	// Normal constructor
			FileReader(const FileReader&) = delete;		// Do not allow use of the copy constructor
			~FileReader();								// Destructor

			// Operators
			FileReader& operator=(const FileReader&) = delete;	// Do not allow use of copy-assignment
			
			// Public methods
			bool IsOpen();	// Checks whether the requested file was opened properly
			MSDFileObject ReadObject();
			ObjectGroup ObjectGroupType() const {return this->objectGroupType;};
			std::string ObjectGroupName() const {return this->objectGroupName;};
			
			// Static methods
			static std::vector<std::string> GetFileList();
			static std::map<std::string, std::string> GetDefinitions();
			static bool AddDefinition(const std::string&, const std::string& _value = "");
	};
}

#endif