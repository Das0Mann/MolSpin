/////////////////////////////////////////////////////////////////////////
// BasicTask implementation (RunSection module)
//
// NOTE: This is the implementation for the base class used for task
// classes. Be careful when changing this class implementation as it will
// affect all task classes!
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ObjectParser.h"
#include "RunSection.h"
#include "BasicTask.h"
#include "Settings.h"
#include "Spin.h"
#include "Interaction.h"
#include "Transition.h"
#include "SpinSystem.h"

namespace RunSection
{
	// -----------------------------------------------------
	// BasicTask Constructors and Destructor
	// -----------------------------------------------------
	BasicTask::BasicTask(const MSDParser::ObjectParser& _properties, const RunSection& _runsection)	: properties( std::make_shared<MSDParser::ObjectParser>(_properties) ), runsection(_runsection), output(),
																										isValid(false), isValidated(false), scalars(nullptr), vectors(nullptr), usedScalars(), usedVectors()
	{
		if(!_runsection.noCalculations)
		{
			bool append;
			std::string str;
			
			// Check whether a logfile is specified for the calculation method
			if(this->properties->Get("logfile", str)
				|| this->properties->Get("log", str)
				|| this->properties->Get("output", str))
			{
				// Check whether we should append to or overwrite previous files
				if(!this->properties->Get("append",append)
					&& !this->properties->Get("appendlog",append))
					append = false;
				
				// The commandline "--append" option
				if(_runsection.overruleAppend) {append = true;}
				
				// Attempt to create the stream
				if(!this->output.SetLogFile(str, append))
				{
					// If it was not possible to create the stream, terminate
					std::cout << "ERROR: Failed to set logfile of task \"" << this->Name() << "\"!" << std::endl;
					std::cout << "----> Requested file \"" << str << "\"." << std::endl;
					std::cout << "----> Make sure that the folder exists and you have permission to write to the file." << std::endl;
					std::exit(1);
				}
			}
			
			// Check whether a data file is specified for the calculation method
			if(this->properties->Get("datafile", str)
				|| this->properties->Get("data", str))
			{
				// Check whether we should append to or overwrite previous files
				if(!this->properties->Get("append",append)
					&& !this->properties->Get("appenddata",append))
					append = false;
				
				// The commandline "--append" option
				if(_runsection.overruleAppend) {append = true;}
				
				// Attempt to create the stream
				if(!this->output.SetDataFile(str, append))
				{
					// If it was not possible to create the stream, terminate
					std::cout << "ERROR: Failed to set datafile of task \"" << this->Name() << "\"!" << std::endl;
					std::cout << "----> Requested file \"" << str << "\"." << std::endl;
					std::cout << "----> Make sure that the folder exists and you have permission to write to the file." << std::endl;
					std::exit(1);
				}
			}
		}
	}
	
	BasicTask::~BasicTask()
	{
	}
	
	// -----------------------------------------------------
	// Does some extra initialization and then calls Validate
	// -----------------------------------------------------
	bool BasicTask::DoValidation()
	{
		std::string str;
		
		// Get the general notification level
		MessageType msgtype = this->runsection.GetSettings()->NotificationLevel();
		
		// If a notification level was specified for this specific task, attempt to parse it
		if(this->properties->Get("notifications", str) || this->properties->Get("notificationlevel", str))
		{
			if(str.compare("details") == 0 || str.compare("showdetails") == 0 || str.compare("all") == 0 || str.compare("verbose") == 0)
			{
				msgtype &= ~NormalMessage;			// Reset normal message type bits
				msgtype |= MessageType_Details;		// Set the bit for the requested type
			}
			else if(str.compare("normal") == 0 || str.compare("standard") == 0)
			{
				msgtype &= ~NormalMessage;			// Reset normal message type bits
				msgtype |= MessageType_Normal;		// Set the bit for the requested type
			}
			else if(str.compare("sparse") == 0 || str.compare("important") == 0 || str.compare("importantonly") == 0)
			{
				msgtype &= ~NormalMessage;			// Reset normal message type bits
				msgtype |= MessageType_Important;	// Set the bit for the requested type
			}
			else if(str.compare("quiet") == 0 || str.compare("silent") == 0 || str.compare("critical") == 0 || str.compare("criticalonly") == 0)
			{
				msgtype &= ~NormalMessage;			// Reset normal message type bits
				msgtype |= MessageType_Critical;	// Set the bit for the requested type
			}
			else
			{
				this->Log(MessageType_Warning) << "Failed to set the notification level! Unrecognized value \"" << str << "\"!" << std::endl;
			}
		}
		
		// The possibility to ignore warnings (not recommended, but possible)
		bool flag;
		if(this->properties->Get("ignorewarnings", flag) || this->properties->Get("nowarn", flag))
		{
			if(flag)
				msgtype &= ~MessageType_Warning;	// Reset the warning flag
			else
				msgtype |= MessageType_Warning;		// Set the warning flag
		}
		
		// The possibility to ignore errors (not recommended, but possible)
		if(this->properties->Get("ignoreerrors", flag) || this->properties->Get("noerr", flag))
		{
			if(flag)
				msgtype &= ~MessageType_Error;	// Reset the error flag
			else
				msgtype |= MessageType_Error;	// Set the error flag
		}
		
		// Set the notification level
		this->output.SetNotificationLevel(msgtype);
		
		// Now call the pure virtual Validate method
		return this->Validate();
	}
	// -----------------------------------------------------
	// BasicTask other private methods
	// -----------------------------------------------------
	// Reset time and trajectory step before running the task
	void BasicTask::ResetTimeAndTrajectoryStep()
	{
		bool trajectoryStepFirst = this->runsection.settings->SetTrajectoryStepBeforeTime();
		for(auto i = this->runsection.systems.cbegin(); i != this->runsection.systems.cend(); i++)
		{
			// Prepare Interaction objects
			auto interactions = (*i)->Interactions();
			for(auto j = interactions.begin(); j != interactions.end(); j++)
			{
				if(trajectoryStepFirst)
				{
					(*j)->SetTrajectoryStep(this->runsection.settings->TrajectoryStep());
					(*j)->SetTime(this->runsection.settings->Time());
				}
				else
				{
					(*j)->SetTime(this->runsection.settings->Time());
					(*j)->SetTrajectoryStep(this->runsection.settings->TrajectoryStep());
				}
			}
		
			// Prepare Transition objects
			auto transitions = (*i)->Transitions();
			for(auto j = transitions.begin(); j != transitions.end(); j++)
			{
				if(trajectoryStepFirst)
				{
					(*j)->SetTrajectoryStep(this->runsection.settings->TrajectoryStep());
					(*j)->SetTime(this->runsection.settings->Time());
				}
				else
				{
					(*j)->SetTime(this->runsection.settings->Time());
					(*j)->SetTrajectoryStep(this->runsection.settings->TrajectoryStep());
				}
			}
		
			// Prepare Spin objects
			auto spins = (*i)->Spins();
			for(auto j = spins.begin(); j != spins.end(); j++)
			{
				if(trajectoryStepFirst)
				{
					(*j)->SetTrajectoryStep(this->runsection.settings->TrajectoryStep());
					(*j)->SetTime(this->runsection.settings->Time());
				}
				else
				{
					(*j)->SetTime(this->runsection.settings->Time());
					(*j)->SetTrajectoryStep(this->runsection.settings->TrajectoryStep());
				}
			}
		}
	}
	// -----------------------------------------------------
	// BasicTask protected methods
	// -----------------------------------------------------
	// This method can be overriden in order to provide a calculation method that runs on supercomputer clusters
	bool BasicTask::RunMPI()
	{
		std::cout << "ERROR: The task \"" << this->Properties()->Name() << "\" does not have an MPI version!" << std::endl;
		return false;
	}
	
	// Method that provides access to the calculation settings
	std::shared_ptr<const Settings> BasicTask::RunSettings() const
	{
		return this->runsection.settings;
	}
	
	// Method that provides access to spin systems
	std::vector<std::shared_ptr<SpinAPI::SpinSystem>> BasicTask::SpinSystems() const	// TODO: Make the shared_ptrs point to const SpinSystems
	{
		return this->runsection.systems;
	}
	
	// Provides access to the properties object for the task
	const std::shared_ptr<MSDParser::ObjectParser>& BasicTask::Properties() const
	{
		return this->properties;
	}
	
	// Returns the log stream
	std::ostream& BasicTask::Log(const MessageType& _msgtype)
	{
		return this->output.Log(_msgtype);
	}
	
	// Returns the data/results stream
	std::ostream& BasicTask::Data()
	{
		return this->output.Data();
	}
	
	// Call the StandardOutput header method
	bool BasicTask::WriteStandardOutputHeader(std::ostream& _stream)
	{
		return this->runsection.WriteOutputHeader(_stream);
	}
	
	// Call the StandardOutput data method
	bool BasicTask::WriteStandardOutput(std::ostream& _stream)
	{
		return this->runsection.WriteOutput(_stream);
	}
	
	// Sets a reference to the requested ActionScalar
	bool BasicTask::Scalar(std::string _name, ActionScalar** _scalar)
	{
		auto i = this->scalars->find(_name);
		if(i != this->scalars->cend())
		{
			// Set the pointer if any
			if(_scalar != nullptr)
			{
				(*_scalar) = const_cast<ActionScalar*>(&(i->second));
				
				// Save the value of the ActionScalar so it can be restored later
				if(this->usedScalars.find(_name) == this->usedScalars.cend())
				{
					this->usedScalars[_name] = std::pair<ActionScalar*, double>(*_scalar, (*_scalar)->Get());
				}
			}
			
			return true;
		}
		
		return false;
	}
	
	// Sets a reference to the requested ActionVector
	bool BasicTask::Vector(std::string _name, ActionVector** _vector)
	{
		auto i = this->vectors->find(_name);
		if(i != this->vectors->cend())
		{
			// Set the pointer if any
			if(_vector != nullptr)
			{
				(*_vector) = const_cast<ActionVector*>(&(i->second));
				
				// Save the value of the ActionVector so it can be restored later
				if(this->usedVectors.find(_name) == this->usedVectors.cend())
				{
					this->usedVectors[_name] = std::pair<ActionVector*, arma::vec>(*_vector, (*_vector)->Get());
				}
			}
			
			return true;
		}
		
		return false;
	}
	// -----------------------------------------------------
	// BasicTask public methods
	// -----------------------------------------------------
	// Runs the task
	bool BasicTask::Run()
	{
		// Make sure the task is valid
		if(!this->IsValid())
		{
			this->Log() << "Warning: Skipping invalid task \"" << this->Name() << "\"." << std::endl;
			return false;
		}
		
		// -----------------------------------------------------
		// Prepare for the task
		// -----------------------------------------------------
		// Set the time and trajectory step on all relevant objects
		this->ResetTimeAndTrajectoryStep();
		
		// Store current value of ActionScalars used by the task
		for(auto i = this->usedScalars.begin(); i != this->usedScalars.end(); i++)
			i->second.second = i->second.first->Get();
		
		// Store current value of ActionVectors used by the task
		for(auto i = this->usedVectors.begin(); i != this->usedVectors.end(); i++)
			i->second.second = i->second.first->Get();
		
		// -----------------------------------------------------
		// Run the task
		// -----------------------------------------------------
		this->Log(MessageType_Important) << "--- Running task \"" << this->Name() << "\" ---" << std::endl;
		bool run_status = this->RunLocal();
		
		// -----------------------------------------------------
		// Clean-up after running the task
		// -----------------------------------------------------
		// Reset ActionScalars to the value stored before the task was run
		for(auto i = this->usedScalars.cbegin(); i != this->usedScalars.cend(); i++)
		{
			if(i->second.first->IsReadonly())
				i->second.first->Reset();
			else
				i->second.first->Set(i->second.second);
		}
		
		// Reset ActionVectors to the value stored before the task was run
		for(auto i = this->usedVectors.cbegin(); i != this->usedVectors.cend(); i++)
		{
			if(i->second.first->IsReadonly())
				i->second.first->Reset();
			else
				i->second.first->Set(i->second.second);
		}
		
		return run_status;
	}
	
	// Checks whether the task is valid
	bool BasicTask::IsValid()
	{
		// Make sure the validation method has been called
		if(!this->isValidated)
		{
			this->isValid = this->DoValidation();
			this->isValidated = true;
			
			// Notify the user if validation fails
			if(!this->isValid)
				this->Log() << "ERROR: Failed to validate task \"" << this->Name() << "\"!" << std::endl;
		}
		
		return this->isValid;
	}
	
	// Returns the name of the task
	std::string BasicTask::Name()
	{
		if(this->properties == nullptr)
			return "undefined";
		
		return this->properties->Name();
	}
	
	// Method to provide BasicTask with access to ActionTargets
	void BasicTask::SetActionTargets(const std::map<std::string, ActionScalar>& _scalars, const std::map<std::string, ActionVector>& _vectors)
	{
		this->scalars = &_scalars;
		this->vectors = &_vectors;
	}
	
	// -----------------------------------------------------
	// Public OutputHandler stream manipulation methods
	// -----------------------------------------------------
	// Change the log stream
	bool BasicTask::SetLogStream(std::ostream& _stream)
	{
		return this->output.SetLogStream(_stream);
	}
	
	// Change the data stream
	bool BasicTask::SetDataStream(std::ostream& _stream)
	{
		return this->output.SetDataStream(_stream);
	}
	// -----------------------------------------------------
}

