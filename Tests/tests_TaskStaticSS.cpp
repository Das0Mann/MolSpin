//////////////////////////////////////////////////////////////////////////////
// MolSpin Unit Testing Module
//
// Tests the Static Superspace method.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#include "TaskStaticSS.h"
//////////////////////////////////////////////////////////////////////////////
// Tests a single calculation without any Actions
bool test_task_staticss_simplemodel()
{
	// Setup objects for the test
	// Spins
	auto spin1 = std::make_shared<SpinAPI::Spin>("electron1", "spin=1/2;tensor=isotropic(2);");
	auto spin2 = std::make_shared<SpinAPI::Spin>("electron2", "spin=1/2;tensor=isotropic(2);");
	auto spin3 = std::make_shared<SpinAPI::Spin>("nucleus1", "spin=1/2;tensor=isotropic(1);");
	auto spin4 = std::make_shared<SpinAPI::Spin>("nucleus2", "spin=1/2;tensor=isotropic(1);");

	// Interactions
	auto interaction1 = std::make_shared<SpinAPI::Interaction>("interaction1", "type=hyperfine;group1=electron1;group2=nucleus1;tensor=isotropic(5e-4);");
	auto interaction2 = std::make_shared<SpinAPI::Interaction>("interaction2", "type=hyperfine;group1=electron2;group2=nucleus2;tensor=anisotropic(1e-4, 1e-4, 1e-3);");
	auto interaction3 = std::make_shared<SpinAPI::Interaction>("interaction3", "type=zeeman;spins=electron1,electron2;field=0 0 5e-5;");

	// States
	auto state1 = std::make_shared<SpinAPI::State>("state1", "spins(electron1,electron2)=|1/2,-1/2>-|-1/2,1/2>;"); // Singlet
	auto state2 = std::make_shared<SpinAPI::State>("state2", "spin(electron1)=|1/2>;spin(electron2)=|1/2>;");	   // |T+>
	auto state3 = std::make_shared<SpinAPI::State>("state3", "");												   // Identity

	// SpinSystem
	auto spinsys = std::make_shared<SpinAPI::SpinSystem>("System");
	spinsys->Add(spin1);
	spinsys->Add(spin2);
	spinsys->Add(spin3);
	spinsys->Add(spin4);
	spinsys->Add(state1);
	spinsys->Add(state2);
	spinsys->Add(state3);
	spinsys->Add(interaction1);
	spinsys->Add(interaction2);
	spinsys->Add(interaction3);
	spinsys->ValidateInteractions();
	std::vector<std::shared_ptr<SpinAPI::SpinSystem>> spinsystems;
	spinsystems.push_back(spinsys);

	// Transition
	auto transition1 = std::make_shared<SpinAPI::Transition>("transition1", "sourcestate=state3;rate=1e-4;", spinsys);
	spinsys->Add(transition1);

	// Add an ObjectParser with settings to the SpinSystem
	auto spinsysParser = std::make_shared<MSDParser::ObjectParser>("spinsyssettings", "initialstate=state1;");
	spinsys->SetProperties(spinsysParser);

	// A RunSection to run the calculation
	RunSection::RunSection rs;
	rs.Add(spinsys);

	// Create a task and get a pointer to it
	std::string taskname = "testtask";
	MSDParser::ObjectParser taskParser(taskname, "type=staticss;");
	rs.Add(MSDParser::ObjectType::Task, taskParser);
	auto task = rs.GetTask(taskname);

	// Set the Log and Data streams to something we can read off within this function
	std::ostringstream logstream;
	std::ostringstream datastream;
	task->SetLogStream(logstream);
	task->SetDataStream(datastream);

	// "Gold values", i.e. correct results to test against.
	std::string value1 = "1 3257.57 1742.56 10000";

	bool isCorrect = true;

	// Perform the test
	isCorrect &= state1->ParseFromSystem(*spinsys);							// Get a valid state object; Singlet
	isCorrect &= state2->ParseFromSystem(*spinsys);							// Get a valid state object; |T+>
	isCorrect &= state3->ParseFromSystem(*spinsys);							// Get a valid state object; Identity
	isCorrect &= ((spinsys->ValidateTransitions(spinsystems)).size() == 0); // Put state object into transition
	isCorrect &= rs.Run(1);													// Run a calculation

	// Remove header from first run
	std::string result_string = datastream.str();
	auto lb = result_string.find("\n");
	if (lb != std::string::npos && lb < result_string.size() - 1)
	{
		result_string.erase(0, lb + 1);
	}

	isCorrect &= equal_doublesfromstring(result_string, value1);

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Like the previous, but with a RotateVector action on the Zeeman field
bool test_task_staticss_simplemodel2()
{
	// Setup objects for the test
	// Spins
	auto spin1 = std::make_shared<SpinAPI::Spin>("electron1", "spin=1/2;tensor=isotropic(2);");
	auto spin2 = std::make_shared<SpinAPI::Spin>("electron2", "spin=1/2;tensor=isotropic(2);");
	auto spin3 = std::make_shared<SpinAPI::Spin>("nucleus1", "spin=1/2;tensor=isotropic(1);");
	auto spin4 = std::make_shared<SpinAPI::Spin>("nucleus2", "spin=1/2;tensor=isotropic(1);");

	// Interactions
	auto interaction1 = std::make_shared<SpinAPI::Interaction>("interaction1", "type=hyperfine;group1=electron1;group2=nucleus1;tensor=isotropic(5e-4);");
	auto interaction2 = std::make_shared<SpinAPI::Interaction>("interaction2", "type=hyperfine;group1=electron2;group2=nucleus2;tensor=anisotropic(1e-4, 1e-4, 1e-3);");
	auto interaction3 = std::make_shared<SpinAPI::Interaction>("interaction3", "type=zeeman;spins=electron1,electron2;field=0 0 5e-5;");

	// States
	auto state1 = std::make_shared<SpinAPI::State>("state1", "spins(electron1,electron2)=|1/2,-1/2>-|-1/2,1/2>;"); // Singlet
	auto state2 = std::make_shared<SpinAPI::State>("state2", "spin(electron1)=|1/2>;spin(electron2)=|1/2>;");	   // |T+>
	auto state3 = std::make_shared<SpinAPI::State>("state3", "");												   // Identity

	// SpinSystem
	auto spinsys = std::make_shared<SpinAPI::SpinSystem>("System");
	spinsys->Add(spin1);
	spinsys->Add(spin2);
	spinsys->Add(spin3);
	spinsys->Add(spin4);
	spinsys->Add(state1);
	spinsys->Add(state2);
	spinsys->Add(state3);
	spinsys->Add(interaction1);
	spinsys->Add(interaction2);
	spinsys->Add(interaction3);
	spinsys->ValidateInteractions();
	std::vector<std::shared_ptr<SpinAPI::SpinSystem>> spinsystems;
	spinsystems.push_back(spinsys);

	// Transition
	auto transition1 = std::make_shared<SpinAPI::Transition>("transition1", "sourcestate=state3;rate=1e-3;", spinsys);
	spinsys->Add(transition1);

	// Add an ObjectParser with settings to the SpinSystem
	auto spinsysParser = std::make_shared<MSDParser::ObjectParser>("spinsyssettings", "initialstate=state1;");
	spinsys->SetProperties(spinsysParser);

	// A RunSection to run the calculation
	RunSection::RunSection rs;
	rs.Add(spinsys);

	// Create a task and get a pointer to it
	std::string taskname = "testtask";
	MSDParser::ObjectParser taskParser(taskname, "type=staticss;");
	rs.Add(MSDParser::ObjectType::Task, taskParser);
	auto task = rs.GetTask(taskname);

	// Create an Action
	MSDParser::ObjectParser actionParser("action1", "type=rotatevector;vector=System.interaction3.field;axis=0 1 0;value=15;");
	rs.Add(MSDParser::ObjectType::Action, actionParser);

	// Set the Log and Data streams to something we can read off within this function
	std::ostringstream logstream;
	std::ostringstream datastream;
	task->SetLogStream(logstream);
	task->SetDataStream(datastream);

	// "Gold values", i.e. correct results to test against. Obtained as in previous test.
	// Obtained by 15 degrees rotation about axis="0 1 0" for the Zeeman field.
	std::string values[25] = {"1 327.212 174.110 1000",
							  "2 321.633 179.687 1000",
							  "3 307.003 194.457 1000",
							  "4 288.503 213.453 1000",
							  "5 271.831 230.932 1000",
							  "6 260.868 242.648 1000",
							  "7 257.143 246.686 1000",
							  "8 260.868 242.648 1000",
							  "9 271.831 230.932 1000",
							  "10 288.503 213.453 1000",
							  "11 307.003 194.457 1000",
							  "12 321.633 179.687 1000",
							  "13 327.212 174.110 1000",
							  "14 321.633 179.687 1000",
							  "15 307.003 194.457 1000",
							  "16 288.503 213.453 1000",
							  "17 271.831 230.932 1000",
							  "18 260.868 242.648 1000",
							  "19 257.143 246.686 1000",
							  "20 260.868 242.648 1000",
							  "21 271.831 230.932 1000",
							  "22 288.503 213.453 1000",
							  "23 307.003 194.457 1000",
							  "24 321.633 179.687 1000",
							  "25 327.212 174.110 1000"};

	bool isCorrect = true;

	// Perform the test
	isCorrect &= state1->ParseFromSystem(*spinsys);							// Get a valid state object; Singlet
	isCorrect &= state2->ParseFromSystem(*spinsys);							// Get a valid state object; |T+>
	isCorrect &= state3->ParseFromSystem(*spinsys);							// Get a valid state object; Identity
	isCorrect &= ((spinsys->ValidateTransitions(spinsystems)).size() == 0); // Put state object into transition
	for (unsigned int i = 0; i < sizeof(values) / sizeof(values[0]); i++)
	{
		isCorrect &= rs.Run(i + 1);

		// Remove header from first run
		std::string result_string = datastream.str();
		auto lb = result_string.find("\n");
		if (lb != std::string::npos && lb < result_string.size() - 1)
		{
			result_string.erase(0, lb + 1);
		}

		isCorrect &= equal_doublesfromstring(result_string, values[i]);
		datastream.str("");
		datastream.clear();
		isCorrect &= rs.Step(i + 2);
	}

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Like the previous, but adds the spins to the spin system in a different
// order, thus testing this feature of the SpinSpace class
// NOTE: This is a test of the SpinSpace::GetState methods, since it is
// identical to the previous test except for the ordering of the spin basis!
bool test_task_staticss_simplemodel2_basisreordering()
{
	// Setup objects for the test
	// Spins
	auto spin1 = std::make_shared<SpinAPI::Spin>("electron1", "spin=1/2;tensor=isotropic(2);");
	auto spin2 = std::make_shared<SpinAPI::Spin>("electron2", "spin=1/2;tensor=isotropic(2);");
	auto spin3 = std::make_shared<SpinAPI::Spin>("nucleus1", "spin=1/2;tensor=isotropic(1);");
	auto spin4 = std::make_shared<SpinAPI::Spin>("nucleus2", "spin=1/2;tensor=isotropic(1);");

	// Interactions
	auto interaction1 = std::make_shared<SpinAPI::Interaction>("interaction1", "type=hyperfine;group1=electron1;group2=nucleus1;tensor=isotropic(5e-4);");
	auto interaction2 = std::make_shared<SpinAPI::Interaction>("interaction2", "type=hyperfine;group1=electron2;group2=nucleus2;tensor=anisotropic(1e-4, 1e-4, 1e-3);");
	auto interaction3 = std::make_shared<SpinAPI::Interaction>("interaction3", "type=zeeman;spins=electron1,electron2;field=0 0 5e-5;");

	// States
	auto state1 = std::make_shared<SpinAPI::State>("state1", "spins(electron1,electron2)=|1/2,-1/2>-|-1/2,1/2>;"); // Singlet
	auto state2 = std::make_shared<SpinAPI::State>("state2", "spin(electron1)=|1/2>;spin(electron2)=|1/2>;");	   // |T+>
	auto state3 = std::make_shared<SpinAPI::State>("state3", "");												   // Identity

	// SpinSystem
	auto spinsys = std::make_shared<SpinAPI::SpinSystem>("System");
	spinsys->Add(spin1);
	spinsys->Add(spin3); // NOTE: The spins are added in a different order
	spinsys->Add(spin2); // This is intended as it tests use of basis reordering in SpinSpace::GetState
	spinsys->Add(spin4);
	spinsys->Add(state1);
	spinsys->Add(state2);
	spinsys->Add(state3);
	spinsys->Add(interaction1);
	spinsys->Add(interaction2);
	spinsys->Add(interaction3);
	spinsys->ValidateInteractions();
	std::vector<std::shared_ptr<SpinAPI::SpinSystem>> spinsystems;
	spinsystems.push_back(spinsys);

	// Transition
	auto transition1 = std::make_shared<SpinAPI::Transition>("transition1", "sourcestate=state3;rate=1e-3;", spinsys);
	spinsys->Add(transition1);

	// Add an ObjectParser with settings to the SpinSystem
	auto spinsysParser = std::make_shared<MSDParser::ObjectParser>("spinsyssettings", "initialstate=state1;");
	spinsys->SetProperties(spinsysParser);

	// A RunSection to run the calculation
	RunSection::RunSection rs;
	rs.Add(spinsys);

	// Create a task and get a pointer to it
	std::string taskname = "testtask";
	MSDParser::ObjectParser taskParser(taskname, "type=staticss;");
	rs.Add(MSDParser::ObjectType::Task, taskParser);
	auto task = rs.GetTask(taskname);

	// Create an Action
	MSDParser::ObjectParser actionParser("action1", "type=rotatevector;vector=System.interaction3.field;axis=0 1 0;value=15;");
	rs.Add(MSDParser::ObjectType::Action, actionParser);

	// Set the Log and Data streams to something we can read off within this function
	std::ostringstream logstream;
	std::ostringstream datastream;
	task->SetLogStream(logstream);
	task->SetDataStream(datastream);

	// "Gold values", i.e. correct results to test against. Obtained as in previous test.
	// Obtained by 15 degrees rotation about axis="0 1 0" for the Zeeman field.
	std::string values[25] = {"1 327.212 174.110 1000",
							  "2 321.633 179.687 1000",
							  "3 307.003 194.457 1000",
							  "4 288.503 213.453 1000",
							  "5 271.831 230.932 1000",
							  "6 260.868 242.648 1000",
							  "7 257.143 246.686 1000",
							  "8 260.868 242.648 1000",
							  "9 271.831 230.932 1000",
							  "10 288.503 213.453 1000",
							  "11 307.003 194.457 1000",
							  "12 321.633 179.687 1000",
							  "13 327.212 174.110 1000",
							  "14 321.633 179.687 1000",
							  "15 307.003 194.457 1000",
							  "16 288.503 213.453 1000",
							  "17 271.831 230.932 1000",
							  "18 260.868 242.648 1000",
							  "19 257.143 246.686 1000",
							  "20 260.868 242.648 1000",
							  "21 271.831 230.932 1000",
							  "22 288.503 213.453 1000",
							  "23 307.003 194.457 1000",
							  "24 321.633 179.687 1000",
							  "25 327.212 174.110 1000"};

	bool isCorrect = true;

	// Perform the test
	isCorrect &= state1->ParseFromSystem(*spinsys);							// Get a valid state object; Singlet
	isCorrect &= state2->ParseFromSystem(*spinsys);							// Get a valid state object; |T+>
	isCorrect &= state3->ParseFromSystem(*spinsys);							// Get a valid state object; Identity
	isCorrect &= ((spinsys->ValidateTransitions(spinsystems)).size() == 0); // Put state object into transition
	for (unsigned int i = 0; i < sizeof(values) / sizeof(values[0]); i++)
	{
		isCorrect &= rs.Run(i + 1);

		// Remove header from first run
		std::string result_string = datastream.str();
		auto lb = result_string.find("\n");
		if (lb != std::string::npos && lb < result_string.size() - 1)
		{
			result_string.erase(0, lb + 1);
		}

		isCorrect &= equal_doublesfromstring(result_string, values[i]);
		datastream.str("");
		datastream.clear();
		isCorrect &= rs.Step(i + 2);
	}

	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Add all the test cases
void AddTaskStaticSSTests(std::vector<test_case> &_cases)
{
	_cases.push_back(test_case("Task StaticSS test 1", test_task_staticss_simplemodel));
	_cases.push_back(test_case("Task StaticSS test 2", test_task_staticss_simplemodel2));
	_cases.push_back(test_case("Task StaticSS test 2 - With spin reordering (tests SpinSpace::GetState for reordering of basis)", test_task_staticss_simplemodel2_basisreordering));
}
//////////////////////////////////////////////////////////////////////////////
