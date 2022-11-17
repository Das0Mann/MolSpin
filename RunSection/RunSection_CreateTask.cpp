/////////////////////////////////////////////////////////////////////////
// CreateTask method of the RunSection class (RunSection module)
// 
// This is where you add new task classes, so that they can be created
// when specified in the input file. Note that this is also where you
// choose the value for the "type" keyword in the input file that will
// be identified with a certain task class.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
// Include task classes here
#include "TaskStaticSS.h"
#include "TaskStaticHSSymmetricDecay.h"
#include "TaskHamiltonianEigenvalues.h"
#include "TaskStaticRPOnlyHSSymDec.h"
#include "TaskStaticSSTimeEvo.h"
#include "TaskDynamicHSTimeEvo.h"
#include "TaskPeriodicSSTimeEvo.h"
#include "TaskPeriodicHSTimeEvo.h"
#include "TaskGammaCompute.h"
#include "TaskMultiStaticSSTimeEvo.h"
#include "TaskMultiDynamicHSTimeEvo.h"

#include "TaskStaticSSRedfield.h"
#include "TaskStaticSSRedfieldSparse.h"
#include "TaskStaticSSRedfieldTimeEvo.h"
#include "TaskStaticSSRedfieldTimeEvoSparse.h"
#include "TaskMultiStaticSSRedfieldTimeEvo.h"
#include "TaskStaticSSSpectra.h"
#include "TaskStaticSSCIDNP.h"
/////////////////////////////////////////////////////////////////////////
namespace RunSection
{
	std::shared_ptr<BasicTask> RunSection::CreateTask(const std::string& _tasktype, const MSDParser::ObjectParser& _obj)
	{
		// The task pointer to be assigned and returned
		std::shared_ptr<BasicTask> task = nullptr;
	
		// Create a task of the proper type
		if(_tasktype.compare("staticss") == 0 || _tasktype.compare("staticivp") == 0) {task = std::make_shared<TaskStaticSS>(_obj, *this);}
		else if(_tasktype.compare("statichs-symmetricdecay") == 0 || _tasktype.compare("statichs") == 0) {task = std::make_shared<TaskStaticHSSymmetricDecay>(_obj, *this);}
		else if(_tasktype.compare("eigenvalues") == 0 || _tasktype.compare("hamiltonianeigenvalues") == 0) {task = std::make_shared<TaskHamiltonianEigenvalues>(_obj, *this);}
		else if(_tasktype.compare("rp-symmetricuncoupled") == 0 || _tasktype.compare("rp-uncoupled") == 0) {task = std::make_shared<TaskStaticRPOnlyHSSymDec>(_obj, *this);}
		else if(_tasktype.compare("staticss-timeevolution") == 0) {task = std::make_shared<TaskStaticSSTimeEvo>(_obj, *this);}
		else if(_tasktype.compare("dynamichs-timeevolution") == 0) {task = std::make_shared<TaskDynamicHSTimeEvo>(_obj, *this);}
		else if(_tasktype.compare("periodicss-timeevolution") == 0) {task = std::make_shared<TaskPeriodicSSTimeEvo>(_obj, *this);}
		else if(_tasktype.compare("periodichs-timeevolution") == 0) {task = std::make_shared<TaskPeriodicHSTimeEvo>(_obj, *this);}
		else if(_tasktype.compare("gamma-compute") == 0) {task = std::make_shared<TaskGammaCompute>(_obj, *this);}
		else if(_tasktype.compare("multistaticss-timeevolution") == 0 || _tasktype.compare("staticss-multisystem") == 0) {task = std::make_shared<TaskMultiStaticSSTimeEvo>(_obj, *this);}
		else if(_tasktype.compare("multidynamichs-timeevolution") == 0 || _tasktype.compare("dynamichs-multisystem") == 0) {task = std::make_shared<TaskMultiDynamicHSTimeEvo>(_obj, *this);}
		
		//NEW (Added by Luca Gerhards): Including BWR-Relaxation Theory Task Class [TaskStaticSSRedfield.cpp] as CreateTask member.
		
		else if(_tasktype.compare("redfield-relaxation") ==0 || _tasktype.compare("Redfield-Relaxation") == 0) {task = std::make_shared<TaskStaticSSRedfield>(_obj, *this);}
		else if(_tasktype.compare("redfield-relaxation-sparse") ==0 || _tasktype.compare("Redfield-Relaxation-Sparse") == 0) {task = std::make_shared<TaskStaticSSRedfieldSparse>(_obj, *this);}
		else if(_tasktype.compare("redfield-relaxation-timeevolution") ==0 || _tasktype.compare("Redfield-Relaxation-Timeevolution") == 0) {task = std::make_shared<TaskStaticSSRedfieldTimeEvo>(_obj, *this);}
		else if(_tasktype.compare("redfield-relaxation-timeevolution-sparse") ==0 || _tasktype.compare("Redfield-Relaxation-Timeevolution-Sparse") == 0) {task = std::make_shared<TaskStaticSSRedfieldTimeEvoSparse>(_obj, *this);}
		else if(_tasktype.compare("multistaticss-redfield-timeevolution") ==0 || _tasktype.compare("staticss-multisystem-redfield-timeevolution") == 0) {task = std::make_shared<TaskMultiStaticSSRedfieldTimeEvo>(_obj, *this);}	
		else if(_tasktype.compare("staticss-spectra") ==0 || _tasktype.compare("StaticSS-Spectra") == 0) {task = std::make_shared<TaskStaticSSSpectra>(_obj, *this);}
		else if(_tasktype.compare("staticss-cidnp") ==0 || _tasktype.compare("StaticSS-CIDNP") == 0) {task = std::make_shared<TaskStaticSSCIDNP>(_obj, *this);}	
		// NOTE: To add a new task class, just add another "else if" here...
		// The string used in the "compare" method is the "type" to be specified in the MolSpin input file
	
		// Return the task instance pointer - or nullptr if the type was not recognized
		return task;
	}
}
/////////////////////////////////////////////////////////////////////////
