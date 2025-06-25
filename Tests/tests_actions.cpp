//////////////////////////////////////////////////////////////////////////////
// MolSpin Unit Testing Module
//
// Unit test functions for the Action classes and ActionTargets.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#include "ActionTarget.h"
#include "ActionAddScalar.h"
#include "ActionMultiplyScalar.h"
#include "ActionAddVector.h"
#include "ActionScaleVector.h"
#include "ActionFibonacciSphere.h"
#include "ActionRotateVector.h"
#include "ActionLogSpace.h"

#include <iostream>
//////////////////////////////////////////////////////////////////////////////
// Tests the ActionTarget alias ActionScalar
// First we define a check-function
bool check_for_test_actiontargets_scalar(const double &_d) { return (_d < 100.0); }
// And then we define the test itself
bool test_actiontargets_scalar()
{
	// Setup objects for the test
	double d = 42.0;
	RunSection::ActionScalar as = RunSection::ActionScalar(d, &check_for_test_actiontargets_scalar);

	bool isCorrect = true;

	// Perform the test
	isCorrect &= equal_double(d, as.Get());
	isCorrect &= as.Set(10.01);
	isCorrect &= equal_double(d, 10.01);
	isCorrect &= equal_double(d, as.Get());
	isCorrect &= !as.Set(100.01); // Check function is false for values >100.0
	isCorrect &= equal_double(d, 10.01);
	isCorrect &= equal_double(d, as.Get());

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the ActionTarget alias ActionVector
// First we define a check-function
bool check_for_test_actiontargets_vector(const arma::vec &_v) { return (!_v.has_nan() && !_v.has_inf() && _v.n_elem == 3); }
// And then we define the test itself
bool test_actiontargets_vector()
{
	// Setup objects for the test
	arma::vec v("1 0 0");
	RunSection::ActionVector av = RunSection::ActionVector(v, &check_for_test_actiontargets_vector);

	bool isCorrect = true;
	arma::vec testvec1("0 42 120");
	arma::vec testvec2("1 1");

	// Perform the test
	isCorrect &= equal_vec(v, av.Get());
	isCorrect &= av.Set(testvec1);
	isCorrect &= equal_vec(v, testvec1);
	isCorrect &= equal_vec(v, av.Get());
	isCorrect &= !av.Set(testvec2); // Check function is false for number of elements != 3
	isCorrect &= equal_vec(v, testvec1);
	isCorrect &= equal_vec(v, av.Get());

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the AddScalar Action
bool test_action_addscalar()
{
	// Setup objects for the test
	double value1 = 42.0;
	double value2 = 62.0;
	double value3 = 82.0;
	double value4 = 102.0;
	double d = value1;
	RunSection::ActionScalar as = RunSection::ActionScalar(d, nullptr); // The ActionScalar (without check function)
	std::map<std::string, RunSection::ActionScalar> asMap;				// ActionScalar map
	std::map<std::string, RunSection::ActionVector> avMap;				// ActionVector map
	asMap.insert(RunSection::NamedActionScalar("testscalar", as));		// Add ActionScalar to the map with the name "testscalar"
	std::string actionname = "test";
	std::string actioncontents = "scalar=testscalar;value=20;";
	MSDParser::ObjectParser parser(actionname, actioncontents);
	RunSection::ActionAddScalar action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// Perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_double(d, value1);
	action.Step(2);
	isCorrect &= equal_double(as.Get(), value2);
	action_ptr->Step(3);
	isCorrect &= equal_double(d, value3);
	action.Step(4);
	isCorrect &= equal_double(as.Get(), value4);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the MultiplyScalar Action
bool test_action_multiplyscalar()
{
	// Setup objects for the test
	double value1 = 42.0;
	double value2 = 63.0;
	double value3 = 94.5;
	double value4 = 141.75;
	double d = value1;
	RunSection::ActionScalar as = RunSection::ActionScalar(d, nullptr); // The ActionScalar (without check function)
	std::map<std::string, RunSection::ActionScalar> asMap;				// ActionScalar map
	std::map<std::string, RunSection::ActionVector> avMap;				// ActionVector map
	asMap.insert(RunSection::NamedActionScalar("testscalar", as));		// Add ActionScalar to the map with the name "testscalar"
	std::string actionname = "test";
	std::string actioncontents = "scalar=testscalar;value=1.5;";
	MSDParser::ObjectParser parser(actionname, actioncontents);
	RunSection::ActionMultiplyScalar action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// Perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_double(d, value1);
	action.Step(2);
	isCorrect &= equal_double(as.Get(), value2);
	action_ptr->Step(3);
	isCorrect &= equal_double(d, value3);
	action.Step(4);
	isCorrect &= equal_double(as.Get(), value4);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the AddVector Action
bool test_action_addvector()
{
	// Setup objects for the test
	arma::vec value1("1 0 0");
	arma::vec value2("1 0 2");
	arma::vec value3("1 0 4");
	arma::vec value4("1 0 6");
	arma::vec v = value1;
	RunSection::ActionVector av = RunSection::ActionVector(v, nullptr); // The ActionVector (without check function)
	std::map<std::string, RunSection::ActionScalar> asMap;				// ActionScalar map
	std::map<std::string, RunSection::ActionVector> avMap;				// ActionVector map
	avMap.insert(RunSection::NamedActionVector("testvec", av));			// Add ActionVector to the map with the name "testvec"
	std::string actionname = "test";
	std::string actioncontents = "vector=testvec;value=2;direction=0 0 100;"; // The "direction" vector should be normalized to 1
	MSDParser::ObjectParser parser(actionname, actioncontents);
	RunSection::ActionAddVector action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// Perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_vec(v, value1);
	action.Step(2);
	isCorrect &= equal_vec(av.Get(), value2);
	action_ptr->Step(3);
	isCorrect &= equal_vec(v, value3);
	action.Step(4);
	isCorrect &= equal_vec(av.Get(), value4);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the ScaleVector Action
bool test_action_scalevector()
{
	// Setup objects for the test
	arma::vec value1("1 0 0");
	arma::vec value2("2 0 0");
	arma::vec value3("4 0 0");
	arma::vec value4("8 0 0");
	arma::vec v = value1;
	RunSection::ActionVector av = RunSection::ActionVector(v, nullptr); // The ActionVector (without check function)
	std::map<std::string, RunSection::ActionScalar> asMap;				// ActionScalar map
	std::map<std::string, RunSection::ActionVector> avMap;				// ActionVector map
	avMap.insert(RunSection::NamedActionVector("testvec", av));			// Add ActionVector to the map with the name "testvec"
	std::string actionname = "test";
	std::string actioncontents = "vector=testvec;value=2;";
	MSDParser::ObjectParser parser(actionname, actioncontents);
	RunSection::ActionScaleVector action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// Perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_vec(v, value1);
	action.Step(2);
	isCorrect &= equal_vec(av.Get(), value2);
	action_ptr->Step(3);
	isCorrect &= equal_vec(v, value3);
	action.Step(4);
	isCorrect &= equal_vec(av.Get(), value4);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the RotateVector Action
bool test_action_rotatevector()
{
	// Setup objects for the test
	arma::vec value1("1 0 0");
	arma::vec value2("0 1 0");
	arma::vec value3("-1 0 0");
	arma::vec value4("0 -1 0");
	arma::vec v = value1;
	RunSection::ActionVector av = RunSection::ActionVector(v, nullptr); // The ActionVector (without check function)
	std::map<std::string, RunSection::ActionScalar> asMap;				// ActionScalar map
	std::map<std::string, RunSection::ActionVector> avMap;				// ActionVector map
	avMap.insert(RunSection::NamedActionVector("testvec", av));			// Add ActionVector to the map with the name "testvec"
	std::string actionname = "test";
	std::string actioncontents = "vector=testvec;value=90;axis=0 0 1;";
	MSDParser::ObjectParser parser(actionname, actioncontents);
	RunSection::ActionRotateVector action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// Perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_vec(v, value1);
	action.Step(2);
	isCorrect &= equal_vec(av.Get(), value2);
	action_ptr->Step(3);
	isCorrect &= equal_vec(v, value3);
	action.Step(4);
	isCorrect &= equal_vec(av.Get(), value4);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the FibonacciSphere Action
bool test_action_fibonaccisphere()
{
	arma::vec value1("0 1 0");
	arma::vec value2("0.707203326 0.4974874372 0.5023641165");		// i = 50;
	arma::vec value3("0.3292513656 -0.005025125628 0.9442289375");	// i = 100
	arma::vec value4("-0.09728003857 -0.9899497487 -0.1026454532"); // i = 198
	arma::vec v = value1;
	RunSection::ActionVector av = RunSection::ActionVector(v, nullptr);
	std::map<std::string, RunSection::ActionScalar> asMap;
	std::map<std::string, RunSection::ActionVector> avMap;
	avMap.insert(RunSection::NamedActionVector("testvec", av));
	std::string actioname = "test";
	std::string actioncontents = "vector=testvec;points=200;";
	MSDParser::ObjectParser parser(actioname, actioncontents);
	RunSection::ActionFibonacciSphere action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_vec(v, value1, 1e-4);
	for (int i = 1; i <= 50; i++)
	{
		action.Step(i);
	}
	isCorrect &= equal_vec(v, value2, 1e-4);
	for (int i = 51; i <= 100; i++)
	{
		action_ptr->Step(i);
	}
	isCorrect &= equal_vec(av.Get(), value3, 1e-4);
	for (int i = 101; i <= 198; i++)
	{
		action.Step(i);
	}
	isCorrect &= equal_vec(v, value4, 1e-4);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Tests the Logspace action
bool test_action_LogSpace()
{
	arma::rowvec logspace = arma::logspace<arma::rowvec>(0, 6, 50);
	double d = logspace[0];
	RunSection::ActionScalar as = RunSection::ActionScalar(d, nullptr); // The ActionScalar (without check function)
	std::map<std::string, RunSection::ActionScalar> asMap;				// ActionScalar map
	std::map<std::string, RunSection::ActionVector> avMap;				// ActionVector map
	asMap.insert(RunSection::NamedActionScalar("testscalar", as));
	std::string actionname = "test";
	std::string actioncontents = "scalar=testscalar;points=50;minvalue=0.0;maxvalue=6.0;";
	MSDParser::ObjectParser parser(actionname, actioncontents);
	RunSection::ActionLogSpace action(parser, asMap, avMap);
	RunSection::Action *action_ptr = &action;

	bool isCorrect = true;

	// perform the test
	isCorrect &= action_ptr->Validate();
	isCorrect &= equal_double(d, logspace[0]);
	for (int i = 1; i <= 10; i++)
	{
		action.Step(i);
	}
	isCorrect &= equal_double(d, logspace[10]);
	for (int i = 11; i <= 20; i++)
	{
		action_ptr->Step(i);
	}
	isCorrect &= equal_double(d, logspace[20]);
	for (int i = 21; i < 50; i++)
	{
		action.Step(i);
	}
	isCorrect &= equal_double(d, logspace[49]);

	// Return the result
	return isCorrect;
}

// Add all the Action classes test cases
void AddActionsTests(std::vector<test_case> &_cases)
{
	_cases.push_back(test_case("RunSection::ActionScalar test with check function", test_actiontargets_scalar));
	_cases.push_back(test_case("RunSection::ActionVector test with check function", test_actiontargets_vector));
	_cases.push_back(test_case("Action AddScalar", test_action_addscalar));
	_cases.push_back(test_case("Action MultiplyScalar", test_action_multiplyscalar));
	_cases.push_back(test_case("Action AddVector", test_action_addvector));
	_cases.push_back(test_case("Action ScaleVector", test_action_scalevector));
	_cases.push_back(test_case("Action RotateVector", test_action_rotatevector));
	_cases.push_back(test_case("Action FibonacciSphere", test_action_fibonaccisphere));
	_cases.push_back(test_case("Action Logspace", test_action_LogSpace));
}
//////////////////////////////////////////////////////////////////////////////
