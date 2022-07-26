//////////////////////////////////////////////////////////////////////////////
// MolSpin Unit Testing Module
// 
// Tests the Static Hilbert Space method with symmetric recombination.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#include "TaskStaticHSSymmetricDecay.h"
//////////////////////////////////////////////////////////////////////////////
// Tests a single calculation without any Actions
bool test_task_statichssymmetricdecay_simplemodel()
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
	auto state1 = std::make_shared<SpinAPI::State>("state1", "spins(electron1,electron2)=|1/2,-1/2>-|-1/2,1/2>;");	// Singlet
	auto state2 = std::make_shared<SpinAPI::State>("state2", "spin(electron1)=|1/2>;spin(electron2)=|1/2>;");		// |T+>
	
	// SpinSystem
	auto spinsys = std::make_shared<SpinAPI::SpinSystem>("System");
	spinsys->Add(spin1);
	spinsys->Add(spin2);
	spinsys->Add(spin3);
	spinsys->Add(spin4);
	spinsys->Add(state1);
	spinsys->Add(state2);
	spinsys->Add(interaction1);
	spinsys->Add(interaction2);
	spinsys->Add(interaction3);
	spinsys->ValidateInteractions();
	
	// Add an ObjectParser with settings to the SpinSystem
	auto spinsysParser = std::make_shared<MSDParser::ObjectParser>("spinsyssettings", "initialstate=state1;");
	spinsys->SetProperties(spinsysParser);
	
	// A RunSection to run the calculation
	RunSection::RunSection rs;
	rs.Add(spinsys);
	
	// Create a task and get a pointer to it
	std::string taskname = "testtask";
	MSDParser::ObjectParser taskParser(taskname, "type=statichs-symmetricdecay;rateconstant=1e-4;");
	rs.Add(MSDParser::ObjectType::Task, taskParser);
	auto task = rs.GetTask(taskname);
	
	// Set the Log and Data streams to something we can read off within this function
	std::ostringstream logstream;
	std::ostringstream datastream;
	task->SetLogStream(logstream);
	task->SetDataStream(datastream);
	
	// "Gold values", i.e. correct results to test against
	// These values were obtained from a different spin dynamics program written before this one,
	// a more simple one that focused only on radical pairs, and which has been thoroughly tested.
	// See e.g.: Scientific Reports volume 6, 36709 (2016), doi:10.1038/srep36709
	// Which used that previous spin dynamics software for the calculations.
	std::string value1 = "1 0.325757 0.174256";
	
	bool isCorrect = true;
	
	// Perform the test
	isCorrect &= state1->ParseFromSystem(*spinsys);	// Get a valid state object; Singlet
	isCorrect &= state2->ParseFromSystem(*spinsys);	// Get a valid state object; |T+>
	isCorrect &= rs.Run(1);							// Run a calculation
	
	// Remove header from first run
	std::string result_string = datastream.str();
	auto lb = result_string.find("\n");
	if(lb != std::string::npos)
	{
		result_string.erase(0,lb+1);
	}
	
	isCorrect &= equal_doublesfromstring(result_string, value1);
	
	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Like the previous, but with a RotateVector action on the Zeeman field
bool test_task_statichssymmetricdecay_simplemodel2()
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
	auto state1 = std::make_shared<SpinAPI::State>("state1", "spins(electron1,electron2)=|1/2,-1/2>-|-1/2,1/2>;");	// Singlet
	auto state2 = std::make_shared<SpinAPI::State>("state2", "spin(electron1)=|1/2>;spin(electron2)=|1/2>;");		// |T+>
	
	// SpinSystem
	auto spinsys = std::make_shared<SpinAPI::SpinSystem>("System");
	spinsys->Add(spin1);
	spinsys->Add(spin2);
	spinsys->Add(spin3);
	spinsys->Add(spin4);
	spinsys->Add(state1);
	spinsys->Add(state2);
	spinsys->Add(interaction1);
	spinsys->Add(interaction2);
	spinsys->Add(interaction3);
	spinsys->ValidateInteractions();
	
	// Add an ObjectParser with settings to the SpinSystem
	auto spinsysParser = std::make_shared<MSDParser::ObjectParser>("spinsyssettings", "initialstate=state1;");
	spinsys->SetProperties(spinsysParser);
	
	// A RunSection to run the calculation
	RunSection::RunSection rs;
	rs.Add(spinsys);
	
	// Create a task and get a pointer to it
	std::string taskname = "testtask";
	MSDParser::ObjectParser taskParser(taskname, "type=statichs-symmetricdecay;rate=1e-3;");
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
	std::string values[25] = {"1 0.327212 0.174110",
							"2 0.321633 0.179687",
							"3 0.307003 0.194457",
							"4 0.288503 0.213453",
							"5 0.271831 0.230932",
							"6 0.260868 0.242648",
							"7 0.257143 0.246686",
							"8 0.260868 0.242648",
							"9 0.271831 0.230932",
							"10 0.288503 0.213453",
							"11 0.307003 0.194457",
							"12 0.321633 0.179687",
							"13 0.327212 0.174110",
							"14 0.321633 0.179687",
							"15 0.307003 0.194457",
							"16 0.288503 0.213453",
							"17 0.271831 0.230932",
							"18 0.260868 0.242648",
							"19 0.257143 0.246686",
							"20 0.260868 0.242648",
							"21 0.271831 0.230932",
							"22 0.288503 0.213453",
							"23 0.307003 0.194457",
							"24 0.321633 0.179687",
							"25 0.327212 0.174110"};
	
	bool isCorrect = true;
	
	// Perform the test
	isCorrect &= state1->ParseFromSystem(*spinsys);	// Get a valid state object; Singlet
	isCorrect &= state2->ParseFromSystem(*spinsys);	// Get a valid state object; |T+>
	for(unsigned int i = 0; i < sizeof(values)/sizeof(values[0]); i++)
	{
		isCorrect &= rs.Run(i+1);
		
		// Remove header from first run
		std::string result_string = datastream.str();
		auto lb = result_string.find("\n");
		if(lb != std::string::npos && lb < result_string.size() - 1)
		{
			result_string.erase(0,lb+1);
		}
		
		isCorrect &= equal_doublesfromstring(result_string, values[i]);
		datastream.str("");
		datastream.clear();
		isCorrect &= rs.Step(i+2);
	}
	
	// Return the result
	return isCorrect;
}
//////////////////////////////////////////////////////////////////////////////
// Add all the test cases
void AddTaskStaticHSSymmetricDecayTests(std::vector<test_case>& _cases)
{
	_cases.push_back(test_case("Task StaticHS-SymmetricDecay test 1", test_task_statichssymmetricdecay_simplemodel));
	_cases.push_back(test_case("Task StaticHS-SymmetricDecay test 2", test_task_statichssymmetricdecay_simplemodel2));
}
//////////////////////////////////////////////////////////////////////////////
