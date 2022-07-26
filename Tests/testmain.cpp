//////////////////////////////////////////////////////////////////////////////
// MolSpin Unit Testing Module
// 
// Note: Strictly speaking the tests here are not unit tests, since some test
// functions contain more than a single test, and since most classes/methods/
// functions are not tested in complete isolation. But the tests are still
// testing various crucial functionalities.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <armadillo>

#include "Tensor.h"
#include "RunSection.h"
//////////////////////////////////////////////////////////////////////////////
using test_ptr = bool (*)();						// Function pointer alias
using test_case = std::pair<std::string,test_ptr>;	// Function pointer and name
//////////////////////////////////////////////////////////////////////////////
// Include various testing functions
#include "assertfunctions.cpp"
//////////////////////////////////////////////////////////////////////////////
// Files with the other test functions
#include "tests_spinapi.cpp"
#include "tests_msdparser.cpp"
#include "tests_actions.cpp"
#include "tests_TaskStaticHSSymmetricDecay.cpp"
#include "tests_TaskStaticSS.cpp"
#include "tests_TaskStaticRPOnlyHSSymDec.cpp"
//////////////////////////////////////////////////////////////////////////////
// A simple test to test the test module itself
bool this_is_a_test_of_the_test_module()
{
	return true;
}
//////////////////////////////////////////////////////////////////////////////
int main(int argc,char** argv)
{
	std::cout << "# -------------------------------------------------------" << std::endl;
	std::cout << "# Molecular Spin Dynamics" << std::endl;
	std::cout << "# " << std::endl;
	std::cout << "# Developed 2017-2019 by Claus Nielsen and 2021-2022 by Luca Gerhards." << std::endl;
	std::cout << "# (c) Quantum Biology and Computational Physics Group," << std::endl;
	std::cout << "# Carl von Ossietzky University of Oldenburg." << std::endl;
	std::cout << "# For more information see www.molspin.eu" << std::endl;
	std::cout << "# -------------------------------------------------------" << std::endl;
	std::cout << "# Unit Testing Module" << std::endl;
	std::cout << "# -------------------------------------------------------" << std::endl;
	
	// Collections of test cases
	std::vector<test_case> cases;
	std::vector<test_case> failed_cases;
	
	// Add test cases to the list
	cases.push_back(test_case("Test module", this_is_a_test_of_the_test_module));
	AddSpinAPITests(cases);
	AddMSDParserTests(cases);
	AddActionsTests(cases);
	AddTaskStaticHSSymmetricDecayTests(cases);
	AddTaskStaticSSTests(cases);
	AddTaskStaticRPOnlyHSSymDecTests(cases);
	
	// Loop through all test cases and test them
	for(auto i = cases.cbegin(); i != cases.cend(); i++)
	{
		std::cout << "Running test \"" << i->first << "\" ......... " << std::flush;
		if(i->second())
		{
			std::cout << "PASSED!";
		}
		else
		{
			std::cout << "FAILED!";
			failed_cases.push_back(*i);
		}
		std::cout << std::endl;
	}
	
	// Get number of passed and failed tests
	int failed = failed_cases.size();
	int passed = cases.size() - failed;
	
	// Get summary at the end
	std::cout << "# -------------------------------------------------------" << std::endl;
	std::cout << "Testing done!\nPassed:  " << passed << "\nFailed:  " << failed << "\nTotal:  " << cases.size() << std::endl;
	
	// Notify which cases has failed
	if(failed > 0)
	{
		std::cout << "\nThere were failed test cases:" << std::endl;
		for(auto i = failed_cases.cbegin(); i != failed_cases.cend(); i++)
			std::cout << " - " << i->first << std::endl;
	}
}
//////////////////////////////////////////////////////////////////////////////
