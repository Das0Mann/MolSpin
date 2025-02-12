/////////////////////////////////////////////////////////////////////////
// Interaction class (SpinAPI Module)
// ------------------
// Base class for interactions.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <memory>
#include <iostream>
#include <Spin.h>
#include "ObjectParser.h"
#include "Interaction.h"

namespace SpinAPI
{
	// -----------------------------------------------------
	// Interaction Constructors and Destructor
	// -----------------------------------------------------
	// The constructor sets up the interaction parameters, but
	// the spin groups are read in the method ParseSpinGroups instead.
	Interaction::Interaction(std::string _name, std::string _contents) : properties(std::make_shared<MSDParser::ObjectParser>(_name, _contents)), couplingTensor(nullptr),
																		 field({0, 0, 0}), dvalue(0.0), evalue(0.0), group1(), group2(), type(InteractionType::Undefined), fieldType(InteractionFieldType::Static), prefactor(1.0), addCommonPrefactor(true), ignoreTensors(false), isValid(true),
																		 trjHasTime(false), trjHasField(false),  trjHasTensor(false), trjHasPrefactor(false), trjTime(0), trjFieldX(0), trjFieldY(0), trjFieldZ(0), trjPrefactor(0),
																		 tdFrequency(1.0), tdPhase(0.0), tdAxis("0 0 1"), tdPerpendicularOscillation(false), tdInitialField({0, 0, 0}),  tensorType(InteractionTensorType::Static), tdTemperature(0.0), tdDamping(0.0), tdRestoring(0.0), tdTimestep(0) ,
																		 tdInitialTensor(3,3, arma::fill::zeros),  tdStdev(0.0), tdMinFreq(0.0), tdMaxFreq(0.0), tdFreqs{}, tdAmps{}, tdPhases{}, tdComponents(0), tdRandOrients(false), tdThetas{}, tdPhis{}, tdCorrTime(0.0), 
																		 tdAmp(0.0), tdDist(false), tdPrintTensor(false), tdPrintField(false),tdSeed(1),tdAutoseed(false), tdGenerator()//, tdFreqs(3, 3, arma::fill::zeros)//, tdFreqs({0,0,0})
	{
		// Is a trajectory specified?
		std::string str;
		if (this->properties->Get("trajectory", str))
		{
			// Get the directory
			std::string directory;
			if (!this->properties->Get("FileReaderDirectory", directory))
			{
				std::cout << "Warning: Failed to obtain input file directory while initializing Interaction " << this->Name() << "!\n";
				std::cout << "         This may lead to problems when trying to localize the trajectory file." << std::endl;
			}

			// Load the trajectory, and check whether something was loaded
			if (this->trajectory.Load(directory + str))
			{
				if (this->trajectory.Length() > 0)
				{
					// Get column numbers
					this->trjHasTime = this->trajectory.HasColumn("time", trjTime);
					this->trjHasPrefactor = this->trajectory.HasColumn("prefactor", trjPrefactor);
					this->trjHasField = this->trajectory.HasColumn("field.x", trjFieldX);
					this->trjHasField &= this->trajectory.HasColumn("field.y", trjFieldY);
					this->trjHasField &= this->trajectory.HasColumn("field.z", trjFieldZ);

				}
			}
			else
			{
				std::cout << "ERROR: Failed to load trajectory \"" << (directory + str) << "\" for Interaction object \"" << this->Name() << "\"!" << std::endl;
			}
		}

		// Get the type of the interaction
		if (this->properties->Get("type", str))
		{

			if (str.compare("onespin") == 0 || str.compare("singlespin") == 0 || str.compare("zeeman") == 0)
			{

				// A onespin interaction should have a field attached, either directly or through a trajectory
				arma::vec inField = arma::zeros<arma::vec>(3);

				this->properties->Get("field", inField);

				if (this->properties->Get("field", inField) || this->trjHasField)
				{
					this->field = inField;
					this->type = InteractionType::SingleSpin;
				}

			}
			else if (str.compare("twospin") == 0 || str.compare("doublespin") == 0 || str.compare("hyperfine") == 0 || str.compare("dipole") == 0)
			{
				this->type = InteractionType::DoubleSpin;
			}
			else if (str.compare("exchange") == 0)
			{
				this->type = InteractionType::Exchange;
			}
			else if (str.compare("zfs") == 0)
			{
				this->type = InteractionType::Zfs;

				double indvalue, inevalue;
				this->Properties()->Get("dvalue", indvalue);
				this->Properties()->Get("evalue", inevalue);

				this->dvalue = indvalue;
				this->evalue = inevalue;
			}
		}

		// If we have a valid interaction type, read the other parameters
		if (this->type != InteractionType::Undefined)
		{	

			Tensor inTensor(0);
			if (this->properties->Get("tensor", inTensor))
			{
				this->couplingTensor = std::make_shared<Tensor>(inTensor);
			}

			this->properties->Get("prefactor", this->prefactor);
			this->properties->Get("commonprefactor", this->addCommonPrefactor);
			this->properties->Get("ignoretensors", this->ignoreTensors);
		}


		// One-spin interactions can have a time-dependent field
		if (this->type == InteractionType::SingleSpin)
		{
			// Do we have a trajectory entry for field? Note that the trajectory may not have a time column, in which case there will be no time dependence
			if (this->trjHasField)
			{
				this->fieldType = InteractionFieldType::Trajectory;
			}

			if (this->properties->Get("fieldtype", str) || this->properties->Get("timedependence", str))
			{

				if (str.compare("broadband") == 0 || str.compare("ougeneral") == 0 ){
					// Random Number Generator Preparation
					std::random_device rand_dev;		// random number generator
					this->tdGenerator.seed(rand_dev()); // random number generator
					this->Properties()->Get("autoseed", this->tdAutoseed);

					if (!this->tdAutoseed)
					{
						std::cout << "Autoseed is off." << std::endl;
						this->Properties()->Get("seed", this->tdSeed);
						if (this->tdSeed != 0)
						{
							this->tdGenerator.seed(this->tdSeed);
							std::cout << "Seed number is " << this->tdSeed << "." << std::endl;
						}
						else
						{
							std::cout << "Undefined seed number! Setting to default of 1." << std::endl;
							std::cout << "# ERROR: undefined seed number! Setting to default of 1." << std::endl;
							this->tdSeed = 1;
						}
					}
				}
				

				if (this->trjHasField)
				{
					std::cout << "Warning: Ignored fieldtype for Interaction \"" << this->Name() << "\"! Using trajectory instead." << std::endl;
				}
				// Figure out what type of time-dependence should be used
				else if (str.compare("linearpolarized") == 0 || str.compare("linearpolarization") == 0 || str.compare("oscillating") == 0)
				{
					this->fieldType = InteractionFieldType::LinearPolarization;
					this->properties->Get("frequency", this->tdFrequency);
					this->properties->Get("phase", this->tdPhase);
					this->properties->Get("printfield", this->tdPrintField);

				}
				else if (str.compare("circularpolarized") == 0 || str.compare("circularpolarization") == 0)
				{
					this->fieldType = InteractionFieldType::CircularPolarization;
					this->properties->Get("frequency", this->tdFrequency);
					this->properties->Get("phase", this->tdPhase);
					this->properties->Get("axis", this->tdAxis);
					this->properties->Get("perpendicularoscillations", this->tdPerpendicularOscillation);
					this->properties->Get("printfield", this->tdPrintField);

				}
				else if (str.compare("broadband") == 0)
				{
					std::cout << "BROADBAND" << std::endl;
					this->fieldType = InteractionFieldType::Broadband;
					this->properties->Get("minfreq", this->tdMinFreq);
					this->properties->Get("maxfreq", this->tdMaxFreq);
					this->properties->Get("stdev", this->tdStdev); //TODO: decide if this stdev should be a scaling ampltiude or actual noise ampltiude
					this->properties->Get("components", this->tdComponents);
					this->properties->Get("randomorientations", this->tdRandOrients);
					this->properties->Get("printfield", this->tdPrintField);
					
					//define all the distributions for broadband noise and fill out vectors 
					std::normal_distribution<double> amp_dist(0.0, this->tdStdev);
					std::uniform_real_distribution<double> phase_dist(0, 2.0 * M_PI);
					std::uniform_real_distribution<double> freq_dist(this->tdMinFreq, this->tdMaxFreq);

					std::vector<double> amps;
					std::vector<double> freqs;
					std::vector<double> phases;

					for(int i_comp=0; i_comp<this->tdComponents; i_comp++){
						double phase = phase_dist(this->tdGenerator);
						double freq = freq_dist(this->tdGenerator);
						double amp = amp_dist(this->tdGenerator);

						phases.push_back(phase);
						freqs.push_back(freq);
						amps.push_back(amp);
					}

					this->tdPhases = phases;
					this->tdFreqs = freqs;
					this->tdAmps = amps;

					if (this->tdRandOrients == true){

						std::vector<double> thetas; 
						std::vector<double> phis;

						std::uniform_real_distribution<double> cos_theta_dist(-1.0, 1.0);
    					std::uniform_real_distribution<double> phi_dist(0.0, 2 * M_PI);

						for(int i_comp=0; i_comp<this->tdComponents; i_comp++){
							double theta_BB = std::acos(cos_theta_dist(this->tdGenerator));
							double phi_BB = phi_dist(this->tdGenerator);

							thetas.push_back(theta_BB);
							phis.push_back(phi_BB);
						}

						this->tdThetas = thetas;
						this->tdPhis = phis;
					}
					//TODO: should have the option to sample random orientations
				}
				else if (str.compare("ougeneral") == 0)
				{
					this->fieldType = InteractionFieldType::OUGeneral;
					this->properties->Get("stdev", this->tdStdev);
					this->properties->Get("correlationtime", this->tdCorrTime);
					this->properties->Get("timestep", this->tdTimestep);
					this->properties->Get("printfield", this->tdPrintField);
				}

				else
				{
					std::cout << "Warning: Unknown fieldtype for Interaction \"" << this->Name() << "\"! Assuming static field." << std::endl;
				}
			}

			// Save the initial field - used only for time-dependent interactions where Actions cannot change the field
			this->tdInitialField = this->field;

			// Set the time to 0 by default
			if (this->HasFieldTimeDependence()){
				this->SetTime(0.0);
			}
				
		}

		// Double-spin interactions can have a time-dependent tensor
		if (this->type == InteractionType::DoubleSpin)
		{

			if (this->properties->Get("tensortype", str))
			{	
				if (str.compare("broadband") == 0 || str.compare("ougeneral") == 0 ){
					// Random Number Generator Preparation
					std::random_device rand_dev;		// random number generator
					this->tdGenerator.seed(rand_dev()); // random number generator
					this->Properties()->Get("autoseed", this->tdAutoseed);

					if (!this->tdAutoseed)
					{
						std::cout << "Autoseed is off." << std::endl;
						this->Properties()->Get("seed", this->tdSeed);
						if (this->tdSeed != 0)
						{
							this->tdGenerator.seed(this->tdSeed);
							std::cout << "Seed number is " << this->tdSeed << "." << std::endl;
						}
						else
						{
							std::cout << "Undefined seed number! Setting to default of 1." << std::endl;
							std::cout << "# ERROR: undefined seed number! Setting to default of 1." << std::endl;
							this->tdSeed = 1;
						}
					}
				}
				
				// Do we have a trajectory entry for tensor? Note that the trajectory may not have a time column, in which case there will be no time dependence
				if (this->trjHasTensor)
				{	
					this->tensorType = InteractionTensorType::Trajectory;
					std::cout << "Warning: Ignored tensortype for Interaction \"" << this->Name() << "\"! Using trajectory instead." << std::endl;
				}

				// Figure out what type of time-dependence should be used
				else if (str.compare("monochromatic") == 0)
				{
					this->tensorType = InteractionTensorType::Monochromatic;
					this->properties->Get("frequency", this->tdFrequency);
					this->properties->Get("phase", this->tdPhase);
					this->properties->Get("amplitude", this->tdAmp);
					this->properties->Get("modulatedistance", this->tdDist);
					this->properties->Get("printtensor", this->tdPrintTensor);
				}
				else if (str.compare("broadband") == 0)
				{
					this->tensorType = InteractionTensorType::Broadband;
					this->properties->Get("minfreq", this->tdMinFreq);
					this->properties->Get("maxfreq", this->tdMaxFreq);
					this->properties->Get("stdev", this->tdStdev);
					this->properties->Get("components", this->tdComponents);
				}

				else if (str.compare("ougeneral") == 0)
				{
					this->tensorType = InteractionTensorType::OUGeneral;
					this->properties->Get("stdev", this->tdStdev);
					this->properties->Get("correlationtime", this->tdCorrTime);
					this->properties->Get("timestep", this->tdTimestep);
					this->properties->Get("modulatedistance", this->tdDist);
					this->properties->Get("printtensor", this->tdPrintTensor);
				}

				// else if (str.compare("ouspring") == 0)
				// {
				// 	this->tensorType = InteractionTensorType::OUSpring;
				// 	this->properties->Get("temperature", this->tdTemperature);
				// 	this->properties->Get("damping", this->tdDamping);
				// 	this->properties->Get("restoring", this->tdRestoring);
				// 	this->properties->Get("timestep", this->tdTimestep);
				// 	this->properties->Get("printtensor", this->tdPrintTensor);

				// }
				
				else
				{
					std::cout << "Warning: Unknown tensortype for Interaction \"" << this->Name() << "\"! Assuming static field." << std::endl;
				}
			
			}

			//Set the time to 0 by default
			if (this->HasTensorTimeDependence())
			{	
				this->tdInitialTensor = couplingTensor->LabFrame();
				this->SetTime(0.0);
			}
		}
	}

	Interaction::Interaction(const Interaction &_interaction) : properties(_interaction.properties), couplingTensor(_interaction.couplingTensor), field(_interaction.field), dvalue(_interaction.dvalue), evalue(_interaction.evalue),
																group1(_interaction.group1), group2(_interaction.group2), type(_interaction.type), fieldType(_interaction.fieldType),
																prefactor(_interaction.prefactor), addCommonPrefactor(_interaction.addCommonPrefactor), ignoreTensors(_interaction.ignoreTensors), isValid(_interaction.isValid),
																trjHasTime(_interaction.trjHasTime), trjHasField(_interaction.trjHasField), trjHasTensor(_interaction.trjHasTensor), trjHasPrefactor(_interaction.trjHasPrefactor), 
																trjTime(_interaction.trjTime), trjFieldX(_interaction.trjFieldX), trjFieldY(_interaction.trjFieldY), trjFieldZ(_interaction.trjFieldZ),
																trjPrefactor(_interaction.trjPrefactor), tdFrequency(_interaction.tdFrequency), tdPhase(_interaction.tdPhase),  tdAxis(_interaction.tdAxis), tdPerpendicularOscillation(_interaction.tdPerpendicularOscillation), 
																tdInitialField(_interaction.tdInitialField),  tensorType(_interaction.tensorType), tdTemperature(_interaction.tdTemperature), 
																tdDamping(_interaction.tdDamping), tdRestoring(_interaction.tdRestoring), tdTimestep(_interaction.tdTimestep), tdInitialTensor(_interaction.tdInitialTensor),
																tdStdev(_interaction.tdStdev), tdMinFreq(_interaction.tdMinFreq), tdMaxFreq(_interaction.tdMaxFreq), tdFreqs(_interaction.tdFreqs), tdAmps(_interaction.tdAmps), tdPhases(_interaction.tdPhases), 
																tdComponents(_interaction.tdComponents), tdRandOrients(_interaction.tdRandOrients), tdThetas(_interaction.tdThetas), tdPhis(_interaction.tdPhis), tdCorrTime(_interaction.tdCorrTime), 
																tdAmp(_interaction.tdAmp), tdDist(_interaction.tdDist), tdPrintTensor(_interaction.tdPrintTensor), tdPrintField(_interaction.tdPrintField),tdSeed(_interaction.tdSeed), tdAutoseed(_interaction.tdAutoseed), tdGenerator(_interaction.tdGenerator)
															
	{
	}

	Interaction::~Interaction()
	{
	}
	// -----------------------------------------------------
	// Operators
	// -----------------------------------------------------
	const Interaction &Interaction::operator=(const Interaction &_interaction)
	{	
		this->properties = std::make_shared<MSDParser::ObjectParser>(*(_interaction.properties));
		this->couplingTensor = _interaction.couplingTensor;
		this->field = _interaction.field;
		this->dvalue = _interaction.dvalue;
		this->evalue = _interaction.evalue;
		this->type = _interaction.type;
		this->fieldType = _interaction.fieldType;
		this->prefactor = _interaction.prefactor;
		this->addCommonPrefactor = _interaction.addCommonPrefactor;
		this->ignoreTensors = _interaction.ignoreTensors;
		this->isValid = _interaction.isValid;

		this->trjHasTime = _interaction.trjHasTime;
		this->trjHasField = _interaction.trjHasField;
		this->trjHasTensor = _interaction.trjHasTensor;
		this->trjHasPrefactor = _interaction.trjHasPrefactor;
		this->trjTime = _interaction.trjTime;
		this->trjFieldX = _interaction.trjFieldX;
		this->trjFieldY = _interaction.trjFieldY;
		this->trjFieldZ = _interaction.trjFieldZ;
		this->trjPrefactor = _interaction.trjPrefactor;

		this->tdFrequency = _interaction.tdFrequency;
		this->tdPhase = _interaction.tdPhase;
		this->tdAxis = _interaction.tdAxis;
		this->tdPerpendicularOscillation = _interaction.tdPerpendicularOscillation;
		this->tdInitialField = _interaction.tdInitialField;
		
		this->tensorType = _interaction.tensorType;
		this->tdTemperature = _interaction.tdTemperature;
		this->tdDamping = _interaction.tdDamping;
		this->tdRestoring = _interaction.tdRestoring;
		this->tdTimestep = _interaction.tdTimestep;
		this->tdInitialTensor = _interaction.tdInitialTensor;

		this->tdStdev = _interaction.tdStdev;
		this->tdMinFreq = _interaction.tdMinFreq;
		this->tdMaxFreq = _interaction.tdMaxFreq;
		this->tdFreqs = _interaction.tdFreqs;
		this->tdAmps = _interaction.tdAmps;
		this->tdPhases = _interaction.tdPhases;
		this->tdComponents = _interaction.tdComponents;
		this->tdRandOrients = _interaction.tdRandOrients;
		this->tdThetas = _interaction.tdThetas;
		this->tdPhis = _interaction.tdPhis;

		this->tdCorrTime = _interaction.tdCorrTime;
		this->tdAmp = _interaction.tdAmp;
		this->tdDist = _interaction.tdDist;
		this->tdPrintTensor = _interaction.tdPrintTensor;
		this->tdPrintField = _interaction.tdPrintField;
		this->tdSeed = _interaction.tdSeed;
		this->tdAutoseed = _interaction.tdAutoseed;
		this->tdGenerator = _interaction.tdGenerator;
		
		return (*this);
	}
	// -----------------------------------------------------
	// Public methods
	// -----------------------------------------------------
	std::string Interaction::Name() const
	{
		return this->properties->Name();
	}

	bool Interaction::IsValid() const
	{
		if(!isValid) { return false; } //by default this function should never return at this point; the only time it should ever return is if the interaction hasn't been assigned a subsystem in a task where subsystems are used 
		
		if (this->type == InteractionType::SingleSpin && !this->group1.empty())
			return true;
		else if (this->type == InteractionType::DoubleSpin && !this->group1.empty() && !this->group2.empty())
			return true;
		else if (this->type == InteractionType::Exchange && !this->group1.empty() && !this->group2.empty())
			return true;
		else if (this->type == InteractionType::Zfs && !this->group1.empty())
			return true;

		return false;
	}
	// -----------------------------------------------------
	// Property methods
	// -----------------------------------------------------
	// Returns the field vector
	const arma::vec Interaction::Field() const
	{
		return this->field;
	}

	// Returns the D value for Zfs
	const double Interaction::Dvalue() const
	{
		return this->dvalue;
	}

	// Returns the E value for Zfs
	const double Interaction::Evalue() const
	{
		return this->evalue;
	}

	// Returns the prefactor value
	const double Interaction::Prefactor() const
	{
		return this->prefactor;
	}

	// Checks whether the interaction is time-dependent
	bool Interaction::HasFieldTimeDependence() const
	{
		if (this->fieldType != InteractionFieldType::Static)
			return true;

		return false;
	}

	bool  Interaction::HasTensorTimeDependence() const
	{
		if (this->tensorType != InteractionTensorType::Static)
			return true;

		return false;
	}


	// Checks whether the interaction is time-dependent
	bool Interaction::HasTimeDependence() const
	{
		if ((this->type == InteractionType::SingleSpin && this->HasFieldTimeDependence()) || (this->trjHasTime && this->trjHasPrefactor) || (this->HasTensorTimeDependence()))
			return true;

		if (this->trjHasTime)
		{
			return true;
		}

		return false;
	}

	//TODO : ADD IN TENSOR TIMEDEPENDENCE CONDITIONAL METHODS
	

	// -----------------------------------------------------
	// Trajectory-related methods
	// -----------------------------------------------------
	// Returns the length of the trajectory, if any
	unsigned int Interaction::TrajectoryLength() const
	{
		return this->trajectory.Length();
	}

	// Set the current state of the Interaction to a specific step in the trajectory (row in the trajectory file)
	bool Interaction::SetTrajectoryStep(unsigned int _step)
	{
		// Make sure that we have a trajectory
		if (this->trajectory.Length() < 1 || _step >= this->trajectory.Length())
		{
			if (this->couplingTensor != nullptr)
				return this->couplingTensor->SetTrajectoryStep(_step); // The tensor might also have a trajectory, so don't return false without asking the tensor...
			else
				return false;
		}

		// Check whether the requested trajectory step is valid
		if (_step < this->trajectory.Length())
		{
			// Get the prefactor at the requested trajectory step
			if (this->trjHasPrefactor)
				this->prefactor = this->trajectory.Get(_step, this->trjPrefactor);

			// Get the field at the requested trajectory step
			if (this->trjHasField && this->fieldType == InteractionFieldType::Trajectory)
			{
				arma::vec tmp = arma::zeros<arma::vec>(3);
				tmp(0) = this->trajectory.Get(_step, this->trjFieldX);
				tmp(1) = this->trajectory.Get(_step, this->trjFieldY);
				tmp(2) = this->trajectory.Get(_step, this->trjFieldZ);
				this->field = tmp;
			}

			// Also set the trajectory step for the tensor, if possible
			// TODO: Consider what should happen if tensor trajectory length is different from interaction object trajectory length
			if (this->couplingTensor != nullptr)
				this->couplingTensor->SetTrajectoryStep(_step);
		}

		return true;
	}

	// Set the current state of the tensor to a specific time defined either in the trajectory (requires a "time" column in the trajectory) or the time-dependence function of the field
	bool Interaction::SetTime(double _time)
	{
		// Set the time for the tensor; the tensor trajectory is independent of the trajectory/time dependence of the interaction object
		if (this->couplingTensor != nullptr &&  this->trjHasTime)
			this->couplingTensor->SetTime(_time);

		// If we have a trajectory, which has a "time" column
		if (this->trajectory.Length() > 0 && this->trjHasTime)
		{
			// Get the row corresponding to the time right after (or at) the requested time
			unsigned int row = 0;
			double timestepAbove = 0.0;
			if (!this->trajectory.FirstRowEqGreaterThan(_time, this->trjTime, row, timestepAbove))
			{
				// If no such time was found, i.e. we are past the last element of the trajectory, use the last step of the trajectory
				this->SetTrajectoryStep(this->trajectory.Length() - 1);
			}
			else if (row == 0)
			{
				// If the specified time is before the first element of the trajectory, use the first step of the trajectory
				this->SetTrajectoryStep((unsigned int)0);
			}
			else
			{
				// Prepare a linear interpolation
				double timestepBelow = this->trajectory.Get(row - 1, trjTime);			  // Get time at previous step
				double l = 1 - (timestepAbove - _time) / (timestepAbove - timestepBelow); // Get linear interpolation factor, i.e. a number between 0 and 1, where 0 = previous time step and 1 = next time step in the trajectory

				// Get a linear interpolation for the prefactor
				if (this->trjHasPrefactor)
					this->prefactor = this->trajectory.Get(row, this->trjPrefactor) * l + this->trajectory.Get(row - 1, this->trjPrefactor) * (1 - l);

				// Get a linear interpolation for the field vector
				if (this->trjHasField && this->fieldType == InteractionFieldType::Trajectory)
				{
					arma::vec tmp = arma::zeros<arma::vec>(3);
					tmp(0) = this->trajectory.Get(row, this->trjFieldX) * l + this->trajectory.Get(row - 1, this->trjFieldX) * (1 - l);
					tmp(1) = this->trajectory.Get(row, this->trjFieldY) * l + this->trajectory.Get(row - 1, this->trjFieldY) * (1 - l);
					tmp(2) = this->trajectory.Get(row, this->trjFieldZ) * l + this->trajectory.Get(row - 1, this->trjFieldZ) * (1 - l);
					this->field = tmp;
				}
			}
		}

		std::string tdFilename = this->properties->Name() + ".mst";
		
		//DO PRINTING OF FIELD AND TENSOR HERE
		if(_time == 0){
			if(this->tdPrintField){
				std::ofstream file;
				file.open(tdFilename);
				file << "time, field.x, field.y, field.z" << std::endl;
				file.close();
			}
		}
		
		// Apply the right time-dependence function if not using a trajectory or static field
		if (this->fieldType == InteractionFieldType::LinearPolarization)
			this->field = FieldTimeDependenceLinearPolarization(this->tdInitialField, _time, this->tdFrequency, this->tdPhase);
		else if (this->fieldType == InteractionFieldType::CircularPolarization)
			this->field = FieldTimeDependenceCircularPolarization(this->tdInitialField, _time, this->tdFrequency, this->tdPhase, this->tdAxis, this->tdPerpendicularOscillation);
		else if (this->fieldType == InteractionFieldType::Broadband)
			this->field = FieldTimeDependenceBroadband(this->tdInitialField, _time, this->tdFreqs, this->tdAmps,  this->tdPhases, this->tdThetas, this->tdPhis, this->tdRandOrients, this->tdComponents);
		else if (this->fieldType == InteractionFieldType::OUGeneral)
			this->field = FieldTimeDependenceOUGeneral(this->tdInitialField, this->field, _time, this->tdTimestep, this->tdStdev,  this->tdCorrTime, this->tdGenerator);
		
		if(this->tdPrintField){
			std::ofstream file;
			file.open(tdFilename, std::ofstream::app);
			file << _time << " " << this->field(0) << " " << this->field(1) << " " << this->field(2) <<  std::endl;
			file.close();
		}

		/////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////// TENSOR TIMEDEP FUNCTIONS IN HERE //////////////////////////////////////

		if(_time == 0){
			//for development purposes
			if(this->tdPrintTensor){
				std::ofstream file;
				file.open(tdFilename);
				file << "time "  << "mat.xx " << "mat.xy " << "mat.xz " << "mat.yx " << "mat.yy " << "mat.yz " << "mat.zx " << "mat.zy " << "mat.zz" << std::endl;
				file.close();
			}
		}
		// if (this->tensorType == InteractionTensorType::SinMat){
		// 	TensorTimeDependenceSinMat(this->tdInitialTensor, _time, this->tdFrequency, this->tdPhase);
		// }
		//TODO: COULD DEFINE RANDOM SEED HERE
		if (this->tensorType == InteractionTensorType::Monochromatic){
			TensorTimeDependenceMonochromatic(this->tdInitialTensor, _time, this->tdFrequency, this->tdPhase, this->tdAmp);
		}
		else if(this->tensorType == InteractionTensorType::Broadband){
			TensorTimeDependenceBroadband(this->tdInitialTensor);
		}
		else if (this->tensorType == InteractionTensorType::OUGeneral){
			TensorTimeDependenceOUGeneral(this->tdInitialTensor, _time, this->tdTimestep, this->tdStdev,  this->tdCorrTime);
		}
		// else if (this->tensorType == InteractionTensorType::OUSpring){
		// 	TensorTimeDependenceOUSpring(this->tdInitialTensor, _time, this->tdTimestep, this->tdTemperature, this->tdDamping, this->tdRestoring, this->tdSeed);
		// }
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		return true;
	}

	// -----------------------------------------------------
	// Access to custom properties
	// -----------------------------------------------------
	std::shared_ptr<const MSDParser::ObjectParser> Interaction::Properties() const
	{
		return this->properties;
	}

	// -----------------------------------------------------
	// Public method to read spins into group1 and group2,
	// returns false if some spins were not found
	// -----------------------------------------------------
	bool Interaction::ParseSpinGroups(const std::vector<spin_ptr> &_spinlist)
	{
		bool createdSpinLists = false;

		if (this->type == InteractionType::SingleSpin || this->type == InteractionType::Zfs)
		{
			// Attempt to get a list of spins from the input file
			std::string str;
			if (!this->properties->Get("spins", str) && !this->properties->Get("group1", str))
				return false;

			createdSpinLists = this->AddSpinList(str, _spinlist, this->group1);
		}
		else if (this->type == InteractionType::DoubleSpin || this->type == InteractionType::Exchange)
		{
			// Attempt to get a list of spins from the input file
			std::string str1;
			std::string str2;
			if (!this->properties->Get("group1", str1) || !this->properties->Get("group2", str2))
				return false;

			createdSpinLists = this->AddSpinList(str1, _spinlist, this->group1);
			createdSpinLists &= this->AddSpinList(str2, _spinlist, this->group2, &(this->group1)); // Cross-check with group1 as a spin can be in only one of the groups
		}

		// Check whether we were successful
		if (!createdSpinLists)
		{
			this->group1.clear();
			this->group2.clear();
			return false;
		}

		return true;
	}

	// Splits the string into spin names, find the corresponding spins in _spinlist and adds them to _group
	// Private helper method for the ParseSpinGroups method
	bool Interaction::AddSpinList(const std::string &_names, const std::vector<spin_ptr> &_spinlist, std::vector<spin_ptr> &_group, const std::vector<spin_ptr> *_crossCheck)
	{
		// Iterate through the spin name (while-loop) with a ',' delimiter
		std::stringstream ss(_names);
		bool spin_was_found;
		std::string str;
		while (std::getline(ss, str, ','))
		{
			spin_was_found = false;

			// Get the spin with the name held by str
			for (auto i = _spinlist.cbegin(); i != _spinlist.cend(); i++)
			{
				if ((*i)->Name().compare(str) == 0)
				{
					// Perform crosscheck if requested
					if (_crossCheck != nullptr && std::find((*_crossCheck).cbegin(), (*_crossCheck).cend(), (*i)) != (*_crossCheck).cend())
						return false;

					// Make sure we don't have duplicates in the spin group
					if (std::find(_group.cbegin(), _group.cend(), (*i)) != _group.cend())
						return false;

					// Add it to the group
					_group.push_back(*i);
					spin_was_found = true;
					break;
				}
			}

			if (!spin_was_found)
				return false;
		}

		return true;
	}
	// -----------------------------------------------------
	// Public spin list methods
	// -----------------------------------------------------
	// Returns a list of spins such that no spins in the list interacts with a spin outside the list (via the current Interaction object)
	std::vector<spin_ptr> Interaction::CompleteSet(const spin_ptr &_spin) const
	{
		std::vector<spin_ptr> result(1); // Only allocate space for a single spin_ptr by default

		// Only the DoubleSpin interaction type couples spins.
		// If the spin is found in either group, return a list containing the set-theoretic union of the groups
		if (this->type == InteractionType::DoubleSpin &&
			(std::find(this->group1.cbegin(), this->group1.cend(), _spin) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), _spin) != this->group2.cend()))
		{
			result.reserve(this->group1.size() + this->group2.size());				   // Reserve space for both groups to avoid more than 1 reallocation
			result.insert(result.begin(), this->group1.cbegin(), this->group1.cend()); // Insert at the beginning of the vector
			result.insert(result.end(), this->group2.cbegin(), this->group2.cend());   // Insert after the previously inserted spins
		}
		else if (this->type == InteractionType::Exchange &&
				 (std::find(this->group1.cbegin(), this->group1.cend(), _spin) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), _spin) != this->group2.cend()))
		{
			result.reserve(this->group1.size() + this->group2.size());				   // Reserve space for both groups to avoid more than 1 reallocation
			result.insert(result.begin(), this->group1.cbegin(), this->group1.cend()); // Insert at the beginning of the vector
			result.insert(result.end(), this->group2.cbegin(), this->group2.cend());   // Insert after the previously inserted spins
		}
		else if (this->type == InteractionType::Zfs &&
				 (std::find(this->group1.cbegin(), this->group1.cend(), _spin) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), _spin) != this->group2.cend()))
		{
			result.reserve(this->group1.size() + this->group2.size());				   // Reserve space for both groups to avoid more than 1 reallocation
			result.insert(result.begin(), this->group1.cbegin(), this->group1.cend()); // Insert at the beginning of the vector
			result.insert(result.end(), this->group2.cbegin(), this->group2.cend());   // Insert after the previously inserted spins
		}
		else
		{
			// If the spin is not in either group, the spin itself forms a complete set with respect to the current interaction object
			result.push_back(_spin);
		}

		return result;
	}

	// Extends the list of spins such that no spins in the list interacts with a spin outside the list (via the current Interaction object)
	// NOTE: The order of the spins in the list is changed if and only if the list is extended.
	bool Interaction::CompleteSet(std::vector<spin_ptr> &_list) const
	{
		bool was_extended = false;

		// Only the DoubleSpin interaction type couples spins
		if (this->type == InteractionType::DoubleSpin || this->type == InteractionType::Exchange)
		{
			bool containsInteractingSpins = false;

			// Loop through all elements in the list
			for (auto i = _list.cbegin(); i != _list.cend(); i++)
			{
				// Check whether any of the spins in the list are coupled to another spin via the current Interaction object
				if (std::find(this->group1.cbegin(), this->group1.cend(), (*i)) != this->group1.cend() || std::find(this->group2.cbegin(), this->group2.cend(), (*i)) != this->group2.cend())
				{
					containsInteractingSpins = true;
					break;
				}
			}

			// There are interacting spins in the list, so we may need to extend it
			if (containsInteractingSpins)
			{
				// First get a temporary list that is the set-theoretic completion of the given list
				std::vector<spin_ptr> tmpvec;											  // Make a temporary list...
				tmpvec.reserve(_list.size() + this->group1.size() + this->group2.size()); // ...and reserve memory to hold the following three inserts
				tmpvec.insert(tmpvec.end(), _list.cbegin(), _list.cend());				  // Insert all elements from the list
				tmpvec.insert(tmpvec.end(), this->group1.cbegin(), this->group1.cend());  // Insert all elements from group1
				tmpvec.insert(tmpvec.end(), this->group2.cbegin(), this->group2.cend());  // Insert all elements from group2
				std::sort(tmpvec.begin(), tmpvec.end());								  // Sort the elements
				auto newEnd = std::unique(tmpvec.begin(), tmpvec.end());				  // Remove duplicates (requires list to be sorted)
				tmpvec.resize(std::distance(tmpvec.begin(), newEnd));					  // Shrink vector size to get rid of removed objects (std::unique cannot do that)

				// Then check if we added new elements in the process
				if (tmpvec.size() > _list.size())
				{
					// Replace the list with the temporary one, if the temporary list was an extension
					_list = tmpvec;
					was_extended = true;
				}
			}
		}

		return was_extended;
	}

	// Public method to check whether a list of spins is complete,
	// i.e. does not interact with any spins outside the list
	bool Interaction::IsComplete(const std::vector<spin_ptr> &_list) const
	{
		unsigned int group1Elements = 0;
		unsigned int group2Elements = 0;

		// Loop through all elements in the list
		for (auto i = _list.cbegin(); i != _list.cend(); i++)
		{
			// Count number of elements from group1 found in the list
			if (std::find(this->group1.cbegin(), this->group1.cend(), (*i)) != this->group1.cend())
				++group1Elements;

			// Count number of elements from group2 found in the list
			if (std::find(this->group2.cbegin(), this->group2.cend(), (*i)) != this->group2.cend())
				++group2Elements;
		}

		// If an element from group1 is present, all elements from group2 must be present, and vice versa
		if (group1Elements > 0 && group2Elements < this->group2.size())
			return false;
		else if (group2Elements > 0 && group1Elements < this->group1.size())
			return false;

		return true;
	}
	// -----------------------------------------------------
	// Methods to create ActionTarget objects
	// -----------------------------------------------------
	// Create ActionVectors
	std::vector<RunSection::NamedActionVector> Interaction::CreateActionVectors(const std::string &_system)
	{
		std::vector<RunSection::NamedActionVector> vectors;

		if (this->IsValid())
		{
			if (this->type == InteractionType::SingleSpin)
			{
				// The field vector
				RunSection::ActionVector fieldVector = RunSection::ActionVector(this->field, &CheckActionVectorInteractionField, this->HasFieldTimeDependence());
				vectors.push_back(RunSection::NamedActionVector(_system + "." + this->Name() + ".field", fieldVector));

				// Allows changing the base field used by time-dependence functions
				if (this->HasFieldTimeDependence() && this->fieldType != InteractionFieldType::Trajectory)
				{
					RunSection::ActionVector initialfieldVector = RunSection::ActionVector(this->tdInitialField, &CheckActionVectorInteractionField);
					vectors.push_back(RunSection::NamedActionVector(_system + "." + this->Name() + ".basefield", initialfieldVector));
				}

				// The axis of circular polarized oscillating fields
				if (this->fieldType == InteractionFieldType::CircularPolarization)
				{
					RunSection::ActionVector oscaxisVector = RunSection::ActionVector(this->tdAxis, &CheckActionVectorInteractionField);
					vectors.push_back(RunSection::NamedActionVector(_system + "." + this->Name() + ".axis", oscaxisVector));
				}
			}
		}

		return vectors;
	}

	// Create ActionScalars
	std::vector<RunSection::NamedActionScalar> Interaction::CreateActionScalars(const std::string &_system)
	{
		std::vector<RunSection::NamedActionScalar> scalars;

		if (this->IsValid())
		{
			// We should always have a scalar for the prefactor
			RunSection::ActionScalar prefactorScalar = RunSection::ActionScalar(this->prefactor, &CheckActionScalarInteractionPrefactor, this->trjHasPrefactor);
			scalars.push_back(RunSection::NamedActionScalar(_system + "." + this->Name() + ".prefactor", prefactorScalar));

			// There may be additional scalars for singlespin interactions
			if (this->type == InteractionType::SingleSpin)
			{
				// Frequency and phase
				if (this->fieldType == InteractionFieldType::LinearPolarization || this->fieldType == InteractionFieldType::CircularPolarization)
				{
					RunSection::ActionScalar frequencyScalar = RunSection::ActionScalar(this->tdFrequency, nullptr);
					scalars.push_back(RunSection::NamedActionScalar(_system + "." + this->Name() + ".frequency", frequencyScalar));

					RunSection::ActionScalar phaseScalar = RunSection::ActionScalar(this->tdPhase, nullptr);
					scalars.push_back(RunSection::NamedActionScalar(_system + "." + this->Name() + ".phase", phaseScalar));
				}
			}
		}

		return scalars;
	}

	//TEST
	void Interaction::TensorTimeDependenceSinMat(arma::mat _m, double _time, double _frequency, double _phase)
	{
		_m(0,0) = _m(0,0);
		_m(0,1) = _m(0,1);
		_m(0,2) = _m(0,2) * cos(_frequency * _time + _phase);
		_m(1,0) = _m(1,0);
		_m(1,1) = _m(1,1);
		_m(1,2) = _m(1,2);
		_m(2,0) = _m(2,0)  * cos(_frequency * _time + _phase);
		_m(2,1) = _m(2,1);
		_m(2,2) = _m(2,2);

		this->couplingTensor->SetTensor(_m);
	}

	//Applies monochromatic noise to each tensor component, scaled by the static value
	void Interaction::TensorTimeDependenceMonochromatic(arma::mat _m, double _time, double _frequency, double _phase, double amp)
	{	
		double A_xx = _m(0,0); double A_yy = _m(1,1); double A_zz = _m(2,2);
		double A_xy = _m(0,1); double A_xz = _m(0,2); double A_yz = _m(1,2);

		if(this->tdDist == true){

			arma::vec eigvals = eig_sym(_m);
			arma::uword max_index = arma::index_max(arma::abs(eigvals));
			double avg_D_factor = eigvals(max_index);

			double avg_dist = std::cbrt(-0.00278/avg_D_factor); //cube root
			double new_dist = avg_dist * (1 + amp * std::sin(_frequency * _time + _phase));

			double new_D_factor = -0.00278 / std::pow(new_dist, 3.0);

			arma::mat tdTensor = _m * (new_D_factor/avg_D_factor);

			A_xx = tdTensor(0,0); A_yy = tdTensor(1,1); A_zz = tdTensor(2,2);
			A_xy = tdTensor(0,1); A_xz = tdTensor(0,2); A_yz = tdTensor(1,2);

			if(this->tdPrintTensor == true){
				std::ofstream file;
				std::string tdFilename = this->properties->Name() + ".mst";
				file.open(tdFilename, std::ofstream::app);
				file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
				file.close();
			}

			this->couplingTensor->SetTensor(tdTensor);

		}
		else{
			A_xx = A_xx * (1 + amp * std::sin(_frequency * _time + _phase));
			A_yy = A_yy * (1 + amp * std::sin(_frequency * _time + _phase));
			A_zz = A_zz * (1 + amp * std::sin(_frequency * _time + _phase));
			A_xy = A_xy * (1 + amp * std::sin(_frequency * _time + _phase));
			A_xz = A_xz * (1 + amp * std::sin(_frequency * _time + _phase));
			A_yz = A_yz * (1 + amp * std::sin(_frequency * _time + _phase));
			
			arma::mat tdTensor = {{A_xx, A_xy, A_xz}, 
									{A_xy, A_yy, A_yz},
									{A_xz, A_yz, A_zz}};

			if(this->tdPrintTensor == true){
				std::ofstream file;
				std::string tdFilename = this->properties->Name() + ".mst";
				file.open(tdFilename, std::ofstream::app);
				file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
				file.close();
			}

			this->couplingTensor->SetTensor(tdTensor);
		}
	}

	void Interaction::TensorTimeDependenceBroadband(arma::mat _m){


	}

	//Applies noise generated using the Ornstein-Uhlenbeck process
	void Interaction::TensorTimeDependenceOUGeneral(arma::mat _m, double _time, double _timestep, double _stdev, double _corrtime)
	{
		arma::mat labTensor = this->couplingTensor->GetTensor();

		std::normal_distribution<double> dist(0.0, 1.0);

		double A_xx = labTensor(0,0); double A_yy = labTensor(1,1); double A_zz = labTensor(2,2);
		double A_xy = labTensor(0,1); double A_xz = labTensor(0,2); double A_yz = labTensor(1,2);

		if(_time == 0){
			
			if(this->tdPrintTensor == true){
				std::ofstream file;
				std::string tdFilename = this->properties->Name() + ".mst";
				file.open(tdFilename, std::ofstream::app);
				file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
				file.close();
			}

			arma::mat tdTensor = {{A_xx, A_xy, A_xz}, 
								{A_xy, A_yy, A_yz},
								{A_xz, A_yz, A_zz}};

			this->couplingTensor->SetTensor(tdTensor);
		}

		else{

			if(this->tdDist == true){//if we vary the distance of a dipolar tensor
				
				//calculate D factor and distance for mean matrix
				arma::vec eigvals_mean = eig_sym(_m);
				arma::uword max_index_mean = arma::index_max(arma::abs(eigvals_mean));
				double mean_D_factor = 3.0 * eigvals_mean(max_index_mean)/4.0;
				double mean_dist = std::abs(std::cbrt(-0.00278/mean_D_factor)); //cube root

				//calculate D factor and distance for the previous timestep
				arma::vec eigvals = eig_sym(labTensor);
				arma::uword max_index = arma::index_max(arma::abs(eigvals));
				double old_D_factor = 3.0 * eigvals(max_index)/4.0;
				double old_dist = std::abs(std::cbrt(-0.00278/old_D_factor)); //cube root

				//generate noise for the new time
				double noise_dist = dist(this->tdGenerator);
				double new_dist = old_dist + _timestep * (mean_dist - old_dist) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_dist;
				double new_D_factor = -0.00278 / std::pow(std::abs(new_dist), 3.0);

				arma::mat tdTensor = labTensor * (new_D_factor/old_D_factor);

				A_xx = tdTensor(0,0); A_yy = tdTensor(1,1); A_zz = tdTensor(2,2);
				A_xy = tdTensor(0,1); A_xz = tdTensor(0,2); A_yz = tdTensor(1,2);

				if(this->tdPrintTensor == true){
					std::ofstream file;
					std::string tdFilename = this->properties->Name() + ".mst";
					file.open(tdFilename, std::ofstream::app);
					file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
					file.close();
					
					std::ofstream file_dist;
					std::string tdFilename_dist = this->properties->Name() + "_Separation.mst";
					file.open(tdFilename_dist, std::ofstream::app);
					file_dist << _time << " " << new_dist << std::endl;
					file_dist.close();
				}

				this->couplingTensor->SetTensor(tdTensor);

			}

			else{

				double noise_xx = dist(this->tdGenerator);
				double noise_yy = dist(this->tdGenerator);
				double noise_zz = dist(this->tdGenerator);
				double noise_xy = dist(this->tdGenerator);
				double noise_xz = dist(this->tdGenerator);
				double noise_yz = dist(this->tdGenerator);

				A_xx = A_xx + _timestep * (_m(0,0) - A_xx) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_xx;
				A_yy = A_yy + _timestep * (_m(1,1) - A_yy) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_yy;
				A_zz = A_zz + _timestep * (_m(2,2) - A_zz) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_zz;
				A_xy = A_xy + _timestep * (_m(0,1) - A_xy) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_xy;
				A_xz = A_xz + _timestep * (_m(0,2) - A_xz) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_xz;
				A_yz = A_yz + _timestep * (_m(1,2) - A_yz) /_corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_yz;

				arma::mat tdTensor = {{A_xx, A_xy, A_xz}, 
										{A_xy, A_yy, A_yz},
										{A_xz, A_yz, A_zz}};

				if(this->tdPrintTensor == true){
					std::ofstream file;
					file.open("Example/standard_examples/OUGeneral.mst", std::ofstream::app);
					file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
					file.close();
				}
				//set the new diagonalised matrix within the Tensor object
				this->couplingTensor->SetTensor(tdTensor);	

			}
					
		}
	}

	// //Applies noise generated using the Ornstein-Uhlenbeck process for a Hookean Spring
	// void Interaction::TensorTimeDependenceOUSpring(arma::mat _m, double _time, double _timestep, double _temperature, double _damping,  double _restoring)
	// {	
	// 	//TODO: check behaviour is as expected
	// 	double k_B = 1.380649e-23;
	// 	double D = (k_B * _temperature) / _damping;
		
	// 	//get the tensor from the previous time point
	// 	arma::mat labTensor = this->couplingTensor->GetTensor();
	
	// 	std::normal_distribution<double> dist(0.0, std::sqrt(2.0 * D * _timestep));

	// 	double A_xx = labTensor(0,0); double A_yy = labTensor(1,1); double A_zz = labTensor(2,2);
	// 	double A_xy = labTensor(0,1); double A_xz = labTensor(0,2); double A_yz = labTensor(1,2);

	// 	double noise_xx = dist(this->tdGenerator);
	// 	double noise_yy = dist(this->tdGenerator);
	// 	double noise_zz = dist(this->tdGenerator);
	// 	double noise_xy = dist(this->tdGenerator);
	// 	double noise_xz = dist(this->tdGenerator);
	// 	double noise_yz = dist(this->tdGenerator);

	// 	//this makes sure the matrix at time=0 is that input in the msd file 
	// 	if(_time == 0){
			
	// 		//for development purposes
	// 		if(this->tdPrintTensor == true){
	// 			std::ofstream file;
	// 			file.open("Example/standard_examples/OUSpring.mst");
	// 			file << "time "  << "mat.xx " << "mat.xy " << "mat.xz " << "mat.yx " << "mat.yy " << "mat.yz " << "mat.zx " << "mat.zy " << "mat.zz" << std::endl;
	// 			file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
	// 			file.close();
	// 		}

	// 		arma::mat tdTensor = {{A_xx, A_xy, A_xz}, 
	// 							{A_xy, A_yy, A_yz},
	// 							{A_xz, A_yz, A_zz}};

	// 		this->couplingTensor->SetTensor(tdTensor);
	// 	}
	// 	else{ //actually generate the noise for times greater than zero

	// 		A_xx = A_xx + (-_restoring * (A_xx - _m(0,0))) * _timestep + noise_xx;
	// 		A_yy = A_yy + (-_restoring * (A_yy - _m(1,1))) * _timestep + noise_yy;
	// 		A_zz = A_zz + (-_restoring * (A_zz - _m(2,2))) * _timestep + noise_zz;
	// 		A_xy = A_xy + (-_restoring * (A_xy - _m(0,1))) * _timestep + noise_xy;
	// 		A_xz = A_xz + (-_restoring * (A_xz - _m(0,2))) * _timestep + noise_xz;
	// 		A_yz = A_yz + (-_restoring * (A_yz - _m(1,2))) * _timestep + noise_yz;

		
	// 		arma::mat tdTensor = {{A_xx, A_xy, A_xz}, 
	// 							{A_xy, A_yy, A_yz},
	// 							{A_xz, A_yz, A_zz}};

	// 		if(this->tdPrintTensor == true){
	// 			std::ofstream file;
	// 			file.open("Example/standard_examples/OUSpring.mst", std::ofstream::app);
	// 			file << _time << " " << A_xx<< " " << A_xy<< " " << A_xz<< " " <<A_xy<< " " << A_yy<< " " << A_yz << " " <<A_xz<< " " << A_yz<< " " << A_zz<< std::endl;
	// 			file.close();
	// 		}
	// 		//set the new diagonalised matrix within the Tensor object
	// 		this->couplingTensor->SetTensor(tdTensor);
	// 	}
	// }

	// Method that calls the methods to generate ActionVectors and ActionScalars and inserts them into the given collections
	void Interaction::GetActionTargets(std::vector<RunSection::NamedActionScalar> &_scalars, std::vector<RunSection::NamedActionVector> &_vectors, const std::string &_system)
	{
		// Get ActionTargets from private methods
		auto scalars = this->CreateActionScalars(_system);
		auto vectors = this->CreateActionVectors(_system);

		// Insert them
		_scalars.insert(_scalars.end(), scalars.begin(), scalars.end());
		_vectors.insert(_vectors.end(), vectors.begin(), vectors.end());

		// If we have a tensor assigned, also add the ActionTargets of the Tensor
		if (this->couplingTensor != nullptr)
			this->couplingTensor->GetActionTargets(_scalars, _vectors, _system + "." + this->Name());
	}


	// -----------------------------------------------------
	// Non-member non-friend methods
	// -----------------------------------------------------
	bool IsValid(const Interaction &_interaction)
	{
		return _interaction.IsValid();
	}

	bool IsStatic(const Interaction &_interaction)
	{
		return !(HasTrajectory(_interaction) | _interaction.HasTimeDependence());
	}

	bool HasTensor(const Interaction &_interaction)
	{
		return (_interaction.CouplingTensor() != nullptr);
	}

	bool HasTrajectory(const Interaction &_interaction)
	{
		if (_interaction.TrajectoryLength() > 0)
			return true;

		return false;
	}

	InteractionType Type(const Interaction &_interaction)
	{
		return _interaction.Type();
	}
	
	// -----------------------------------------------------
	// Non-member non-friend time-dependent field functions
	// -----------------------------------------------------
	// Linear polarized oscillating magnetic field
	arma::vec FieldTimeDependenceLinearPolarization(const arma::vec &_v, double _time, double _frequency, double _phase)
	{
		return (_v * cos(_frequency * _time + _phase));
	}

	// Circularly polarized oscillating magnetic field
	arma::vec FieldTimeDependenceCircularPolarization(const arma::vec &_v, double _time, double _frequency, double _phase, const arma::vec &_axis, bool _perpendicularOscillations)
	{
		// Get the angle for the rotation matrix
		double angle = _time * _frequency + _phase;

		// Helper matrix
		arma::mat W(3, 3, arma::fill::zeros);
		W(0, 1) = -_axis(2);
		W(0, 2) = _axis(1);
		W(1, 2) = -_axis(0);
		W(1, 0) = _axis(2);
		W(2, 0) = -_axis(1);
		W(2, 1) = _axis(0);

		// The Rodrigues rotation matrix
		arma::mat R(3, 3, arma::fill::eye);
		R += sin(angle) * W + (2 * sin(angle / 2.0) * sin(angle / 2.0) * W * W);

		if (_perpendicularOscillations)
		{
			// Get the projection of the field vector in the plane perpendicular to the axis, and apply the rotation
			return R * (_v - dot(_axis, _v) * _axis);
		}

		// Just apply the rotation matrix to the field without projecting onto perpendicular plane
		return R * _v;
	}

	//Broadband noise made up of frequency and components sampled from a uniform distribution, and amplitudes sampled from a normal distribution
	arma::vec FieldTimeDependenceBroadband(const arma::vec &_v, double _time, std::vector<double> _freqs, std::vector<double> _amps,  std::vector<double> _phases, std::vector<double> _thetas, std::vector<double> _phis, bool _randorients,  int _comps)
	{
		
		arma::vec new_field = _v;
	
		//if we include random orientations - each component of the time-dependent field has independent broadband noise. If not,
		//a scalar function applies noise to the field terms defined in the interaction
		if(_randorients == true){

			for(int comp=0; comp<_comps; comp++){ 
				double phase = _phases.at(comp);
				double freq = _freqs.at(comp);
				double amp = _amps.at(comp);
				double theta = _thetas.at(comp);
				double phi = _phis.at(comp);

				arma::vec orient = {std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};

				new_field += orient * amp/sqrt(_comps) * std::sin(freq * _time + phase);
			}	
		}
		else{
			double Asinwdt = 0.0;

			if(_comps > 0){
				for(int comp=0; comp<_comps; comp++){ 
					double phase = _phases.at(comp);
					double freq = _freqs.at(comp);
					double amp = _amps.at(comp);

					Asinwdt += amp/sqrt(_comps) * std::sin(freq * _time + phase);
				}
			}
			else{
				Asinwdt = 1.0;
			}	
			new_field += Asinwdt * _v;
		}

		return new_field;
	}

	arma::vec FieldTimeDependenceOUGeneral(const arma::vec &_v_init, const arma::vec &_v_prev, double _time, double _timestep, double _stdev,  double _corrtime, std::mt19937 _generator){

		std::normal_distribution<double> dist(0.0, 1.0);
		double noise_x = dist(_generator);
		double noise_y = dist(_generator);
		double noise_z = dist(_generator);

		double field_x = _v_prev(0) + _timestep * (_v_init(0) - _v_prev(0)) / _corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_x; 
		double field_y = _v_prev(1) + _timestep * (_v_init(1) - _v_prev(1)) / _corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_y; 
		double field_z = _v_prev(2) + _timestep * (_v_init(2) - _v_prev(2)) / _corrtime + _stdev * std::sqrt(2.0 * _timestep/_corrtime) * noise_z; 
		
		arma::vec newField = {field_x, field_y, field_z};

		return newField;
	}

	
	// -----------------------------------------------------
	// Non-member non-friend ActionTarget Check functions
	// -----------------------------------------------------
	// Checks that a vector set by an Action is valid
	bool CheckActionVectorInteractionField(const arma::vec &_v)
	{
		// A field vector must have 3 components
		if (_v.n_elem != 3)
			return false;

		// Make sure we don't have invalid values
		if (_v.has_nan() || _v.has_inf())
			return false;

		return true;
	}

	// Make sure that the prefactor has a valid value (not NaN or infinite)
	bool CheckActionScalarInteractionPrefactor(const double &_d)
	{
		return std::isfinite(_d);
	}
	// -----------------------------------------------------
}