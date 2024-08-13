/////////////////////////////////////////////////////////////////////////
// Operator class (SpinAPI Module)
// ------------------
// Special operators to be used in some task types.
// 
// Molecular Spin Dynamics Software - developed by Luca Gerhards.
// (c) 2024 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_Pulse
#define MOD_SpinAPI_Pulse

#include <memory>
#include <vector>
#include "SpinAPIDefines.h"
#include "MSDParserfwd.h"
#include "SpinAPIfwd.h"
#include "ActionTarget.h"

namespace SpinAPI
{
	class Pulse
	{
		private:
			// Implementation details
			std::shared_ptr<MSDParser::ObjectParser> properties;
			PulseType type;                                         // Different Pulsetype options are available
            std::vector<spin_ptr> group;							// Spins that are affected by the pulse
			double timestep;                                        // Timestep with which propagation should be done 
            arma::vec rotationaxis;                                  // Direction around which a pulse should rotate
            double angle;                                           // Rotation around a certain direction with angler (deg)
            double pulsetime;                                       // Time of the applied pulse (ns)
            arma::vec field;                                       // Magnetic field strength of pulse (T)
			double frequency;                                        // Frequency of the applied field in HHz;
			double prefactor;										// An optional additional prefactor that can be specified in input file (default = 1.0)
			bool addCommonPrefactor;								// Whether or not to multiply by "g mu_B" for electronic spins, or the equivalent for nuclear spins
			std::vector<double> prefactor_list;							   // Spins that are affected by the pulse

            // Helper method called by ParseSpinGroups
			bool AddSpinList(const std::string&, const std::vector<spin_ptr>&, std::vector<spin_ptr>&, const std::vector<spin_ptr>* _crossCheck = nullptr);

            // Private methods to create ActionTargets
			std::vector<RunSection::NamedActionVector> CreateActionVectors(const std::string&);
			std::vector<RunSection::NamedActionScalar> CreateActionScalars(const std::string&);
			
		public:
			// Constructors / Destructors
			Pulse(std::string, std::string);			// Normal constructor
			Pulse(const Pulse&);						// Copy-constructor
			~Pulse();									// Destructor
			
			// Operator
			const Pulse& operator=(const Pulse&);		// Copy-assignment
			
            // Read spins into group1 and group2, returns false if some spins were not found
			bool ParseSpinGroups(const std::vector<spin_ptr>&);

			// Name and validation
			std::string Name() const;
			bool Validate(const std::vector<std::shared_ptr<SpinAPI::SpinSystem>>&);
			bool IsValid() const;
			
			// Public property methods
			PulseType Type() const {return this->type;};
			std::vector<spin_ptr> Group() const {return this->group;};
			const bool AddCommonPrefactor() const {return this->addCommonPrefactor;};
			const double Timestep() const;
			const double Prefactor() const;
            const arma::vec Rotationaxis() const;
            const double Angle() const;
            const double Pulsetime() const;
            const arma::vec Field() const;
			const double Frequency() const;
			const arma::vec Prefactor_list() const ;
			
			
			// Allow access to custom properties to be used for custom tasks
			std::shared_ptr<const MSDParser::ObjectParser> Properties() const;

            	// Public method for creating ActionTargets
			void GetActionTargets(std::vector<RunSection::NamedActionScalar>&, std::vector<RunSection::NamedActionVector>&, const std::string&);
	};
	
	// Define alias for Operator-pointers
	using pulse_ptr = std::shared_ptr<Pulse>;

	PulseType Type(const Pulse&);
	
	// Non-member non-friend functions
	bool IsValid(const Pulse&);

	// Non-member non-friend functions for ActionTarget validation
	bool CheckActionVectorPulseField(const arma::vec &);
	bool CheckActionScalarPulseScalar(const double &);
}

#endif
