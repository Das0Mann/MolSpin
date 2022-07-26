/////////////////////////////////////////////////////////////////////////
// Object Parser implementations, provides easy access to properties for
// the loaded objects.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "Tensor.h"
#include "ObjectParser.h"
#include "SpinAPIfwd.h"

namespace MSDParser
{
	// -----------------------------------------------------
	// ObjectParser Constructors and Destructor
	// -----------------------------------------------------
	// Initialize the ObjectParser with a string
	ObjectParser::ObjectParser(std::string _name, std::string _contents) : fields(), name(_name)
	{
		int parantheses = 0;
		std::string key("");
		std::string value("");
		bool isKey = true;
		auto i = _contents.cbegin();

		while(i != _contents.cend())
		{
			if(isKey)
			{
				if((*i) == '=')			// Finding '=' means we got the key and then need the value
					isKey = false;
				else if((*i) == ';')	// Finding ';' in a key is wrong, so skip the key
					key.clear();
				else
					key += (*i);
			}
			else
			{
				// Track parantheses such that a value of a matrix string, e.g. "1 2;3 4", also can be read (specifed by "matrix(...)" for Tensor objects)
				if((*i) == '(') {++parantheses;}
				else if((*i) == ')') {--parantheses;}
				
				if((*i) == ';' && parantheses < 1)
				{
					// Make sure that neither the key nor the value is empty before using them
					if(key.size() > 0 && value.size() > 0)
						this->fields[key] = value;
					
					// Clear buffers and reset
					key.clear();
					value.clear();
					isKey = true;
				}
				else
				{
					value += (*i);
				}
			}
			
			i++;
		}
	}
	
	// Copy-constructor
	// Copies the values through the initializer list
	ObjectParser::ObjectParser(const ObjectParser& _objParser) : fields(_objParser.fields), name(_objParser.name)
	{
	}
	
	// Destructor
	ObjectParser::~ObjectParser()
	{
	}
	
	// -----------------------------------------------------
	// ObjectParser public Get method overloads
	// -----------------------------------------------------
	// Attemp to find a keyword matching the given name
	bool ObjectParser::Get(const std::string& _str, std::string& _out) const
	{
		auto i = this->fields.find(_str);
		
		if(i != this->fields.end())
			_out = i->second;
		else
 			return false;
		
		return true;
	}
	
	// Attemp to find a double with the given name
	bool ObjectParser::Get(const std::string& _str, double& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		double tmp;
		
		try
		{
			tmp = std::stod(str);
		}
		catch(const std::exception&) {return false;}
		//vector<double>
		_out = tmp;
		return true;
	}
	
	// Attemp to find an unsigned integer with the given name
	bool ObjectParser::Get(const std::string& _str, unsigned int& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		unsigned int tmp;
		
		try
		{
			if(std::stoi(str) < 0)
				return false;
			tmp = static_cast<unsigned int>(std::stoi(str));
		}
		catch(const std::exception&) {return false;}
		
		_out = tmp;
		return true;
	}
	
	// Attemp to find an integer with the given name
	bool ObjectParser::Get(const std::string& _str, int& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		int tmp;
		
		try
		{
			tmp = std::stoi(str);
		}
		catch(const std::exception&) {return false;}
		
		_out = tmp;
		return true;
	}
	
	// Attemp to find a bool with the given name
	bool ObjectParser::Get(const std::string& _str, bool& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		if(str.compare("true") == 0 || str.compare("yes") == 0 || str.compare("1") == 0)
			_out = true;
		else if(str.compare("false") == 0 || str.compare("no") == 0 || str.compare("0") == 0)
			_out = false;
		else
			return false;
		
		return true;
	}
	
	// Attemp to find a tensor with the given name
	bool ObjectParser::Get(const std::string& _str, SpinAPI::Tensor& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		std::string directory = "";
		this->Get("FileReaderDirectory", directory);
		
		_out = SpinAPI::Tensor(str, directory);
		return true;
	}
	
	// Attemp to find a vector with the given name
	bool ObjectParser::Get(const std::string& _str, arma::vec& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		// Parse string and make sure that we have a 3D vector
		arma::vec tmpVec = arma::vec(str);
		if(tmpVec.n_elem != 3)
			return false;
		
		_out = tmpVec;
		return true;
	}
	
	// -----------------------------------------------------
	// ObjectParser public GetList method overloads
	// - Similar to Get, but allows multiple values to be specifed
	// -----------------------------------------------------
	// Attemp to find a keyword matching the given name
	bool ObjectParser::GetList(const std::string& _str, std::vector<std::string>& _out, const char& _delimiter) const
	{
		auto i = this->fields.find(_str);
		
		// If the keyword was found
		if(i != this->fields.end())
		{
			
			// Get all the values separated by the delimiter
			std::istringstream stream(i->second);
			for(std::string s; std::getline(stream, s, _delimiter); )
				_out.push_back(s);
			
			return true;
		}
		
		return false;
	}
	
	// Attemp to find a list of doubles with the given name
	bool ObjectParser::GetList(const std::string& _str, std::vector<double>& _out, const char& _delimiter) const
	{
		// Get a list of strings from the keyword
		std::vector<std::string> strs;
		if(!this->GetList(_str, strs, _delimiter))
			return false;
		
		std::vector<double> tmp;
		tmp.reserve(strs.size());
		
		try
		{
			// Attemp to parse the values
			for(auto i = strs.cbegin(); i != strs.cend(); i++)
				tmp.push_back( std::stod(*i) );
		}
		catch(const std::exception&) {return false;}
		
		_out = tmp;
		return true;
	}
	
	// Attemp to find a list of unsigned integers with the given keyword
	bool ObjectParser::GetList(const std::string& _str, std::vector<unsigned int>& _out, const char& _delimiter) const
	{
		// Get a list of strings from the keyword
		std::vector<std::string> strs;
		if(!this->GetList(_str, strs, _delimiter))
			return false;
		
		std::vector<unsigned int> tmp;
		tmp.reserve(strs.size());
		
		try
		{
			// Attemp to parse the values
			for(auto i = strs.cbegin(); i != strs.cend(); i++)
			{
				if(std::stoi(*i) < 0)
					return false;
				tmp.push_back( static_cast<unsigned int>(std::stoi(*i)) );
			}
		}
		catch(const std::exception&) {return false;}
		
		_out = tmp;
		return true;
	}
	
	// Attemp to find a list of integers with the given keyword
	bool ObjectParser::GetList(const std::string& _str, std::vector<int>& _out, const char& _delimiter) const
	{
		// Get a list of strings from the keyword
		std::vector<std::string> strs;
		if(!this->GetList(_str, strs, _delimiter))
			return false;
		
		std::vector<int> tmp;
		tmp.reserve(strs.size());
		
		try
		{
			// Attemp to parse the values
			for(auto i = strs.cbegin(); i != strs.cend(); i++)
				tmp.push_back( std::stoi(*i) );
		}
		catch(const std::exception&) {return false;}
		
		_out = tmp;
		return true;
	}
	
	// Attemp to find a list of vectors with the given keyword
	bool ObjectParser::GetList(const std::string& _str, std::vector<arma::vec>& _out, const char& _delimiter) const
	{
		// Get a list of strings from the keyword
		std::vector<std::string> strs;
		if(!this->GetList(_str, strs, _delimiter))
			return false;
		
		std::vector<arma::vec> tmp;
		tmp.reserve(strs.size());
		
		// Parse string and make sure that we have a 3D vector
		for(auto i = strs.cbegin(); i != strs.cend(); i++)
		{
			arma::vec tmpVec = arma::vec(*i);
			if(tmpVec.n_elem != 3)
				return false;
			tmp.push_back(tmpVec);
		}
		
		_out = tmp;
		return true;
	}
	
	// -----------------------------------------------------
	// Specialized Get method to get a spin quantum number
	// -----------------------------------------------------
	// Attempt to find a spin quantum number with the given name
	bool ObjectParser::GetSpin(const std::string& _str, unsigned int& _out) const
	{
		std::string str("");
		if(!this->Get(_str,str))
			return false;
		
		unsigned int tmp;
		
		try
		{
			int tmp2 = -1;
			
			// Check whether the spin is written in the form "1/2", "2/2", "3/2", etc.
			if(str.size() > 2 && (*(str.cend()-1) == '2' && *(str.cend()-2) == '/'))
			{
				tmp2 = std::stoi(str.substr(0,str.size()-2).c_str());
			}
			else
			{
				// Use a factor of two if the spin quantum number is not specified with "/2"
				tmp2 = 2 * std::stoi(str.c_str());
			}
			
			
			if(std::stoi(str) < 0)
				return false;
			tmp = static_cast<unsigned int>(tmp2);
		}
		catch(const std::exception&) {return false;}
		
		_out = tmp;
		return true;
	}
	
	// -----------------------------------------------------
	// ObjectParser public GetFunction method
	// -----------------------------------------------------
	std::vector<std::pair<std::string, std::string>> ObjectParser::GetFunction(std::string _keyword) const
	{
		std::vector<std::pair<std::string, std::string>> result;
		
		// Loop through all the fields
		for(auto i = this->fields.cbegin(); i != this->fields.cend(); i++)
		{
			// Search for parantheses
			auto pstart = i->first.find_first_of("(");
			if(pstart != std::string::npos)
			{
				// If the keyword was found, put it into the results collection
				if(i->first.substr(0,pstart).compare(_keyword) == 0)
				{
					auto pend = i->first.find_first_of(")");
					result.push_back(std::pair<std::string, std::string>(i->first.substr(pstart+1,pend-pstart-1), i->second));
				}
			}
		}
		
		return result;
	}
	// -----------------------------------------------------
	// ObjectParser other public methods
	// -----------------------------------------------------
	std::string ObjectParser::Name() const
	{
		return this->name;
	}
	// -----------------------------------------------------
}

