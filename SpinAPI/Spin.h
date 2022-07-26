/////////////////////////////////////////////////////////////////////////
// Spin class (SpinAPI Module)
// ------------------
// Basic spin class to represent e.g. a radical.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_Spin
#define MOD_SpinAPI_Spin

#include <map>
#include <memory>
#include <armadillo>
#include "Tensor.h"
#include "SpinAPIDefines.h"
#include "MSDParserfwd.h"
#include "ActionTarget.h"

namespace SpinAPI
{
	class Spin
	{
		private:
			// Implementation
			Tensor tensor;
			SpinType type;
			std::shared_ptr<MSDParser::ObjectParser> properties;	// Use a pointer to the object to minimize compilation dependencies
			unsigned int s;											// Spin quantum number in units of 1/2*hbar, i.e. 1 for an electron or proton
			arma::vec quantizationAxis1;	// The quantization axes of the spin
			arma::vec quantizationAxis2;
			arma::vec quantizationAxis3;
			Trajectory trajectory;

			bool trjHasTime;
			bool trjHasQAxis1;
			bool trjHasQAxis2;
			bool trjHasQAxis3;
			unsigned int trjTime;
			unsigned int trjQAxis1X;
			unsigned int trjQAxis1Y;
			unsigned int trjQAxis1Z;
			unsigned int trjQAxis2X;
			unsigned int trjQAxis2Y;
			unsigned int trjQAxis2Z;
			unsigned int trjQAxis3X;
			unsigned int trjQAxis3Y;
			unsigned int trjQAxis3Z;

			// Private methods
			void OrthonormalizeQAxes();					// Makes axes orthonormal
			
			// Collections of spin matrices (such that they are only created once)
			static std::map<unsigned int,std::shared_ptr<arma::sp_cx_mat>> SxCollection;
			static std::map<unsigned int,std::shared_ptr<arma::sp_cx_mat>> SyCollection;
			static std::map<unsigned int,std::shared_ptr<arma::sp_cx_mat>> SzCollection;
			static std::map<unsigned int,std::shared_ptr<arma::sp_cx_mat>> SpCollection;
			static std::map<unsigned int,std::shared_ptr<arma::sp_cx_mat>> SmCollection;
			
			// Methods to get the spin matrices
			static std::shared_ptr<arma::sp_cx_mat> SxFromCollection(unsigned int);
			static std::shared_ptr<arma::sp_cx_mat> SyFromCollection(unsigned int);
			static std::shared_ptr<arma::sp_cx_mat> SzFromCollection(unsigned int);
			static std::shared_ptr<arma::sp_cx_mat> SpFromCollection(unsigned int);
			static std::shared_ptr<arma::sp_cx_mat> SmFromCollection(unsigned int);
		
		public:
			// Constructors / Destructors
			Spin(std::string, std::string);	// Normal constructor
			Spin(const Spin&);				// Copy-constructor
			~Spin();						// Destructor
			
			// Operators
			const Spin& operator=(const Spin&);	// Copy-assignment
			
			// Name and validation
			std::string Name() const;
			bool IsValid();
			
			// Public methods
			SpinType Type() const {return this->type;}
			int Multiplicity() const {return (this->s + 1);};		// Spin multiplicity (2s + 1)
			int S() const {return static_cast<int>(this->s);};		// Spin quantum number in units of 1/2*hbar
			const Tensor& GetTensor() const {return this->tensor;};	// Get G-tensor
			bool SetTrajectoryStep(unsigned int);					// Updates g-tensor from its trajectory (if any)
			bool SetTime(double);									// Updates g-tensor from its trajectory (if any)
			unsigned int TrajectoryLength() const;					// Number of steps in the trajectory for g-tensor (0 means no trajectory)
			arma::vec QuantizationAxis1() const;					// First quantization axis of the spin
			arma::vec QuantizationAxis2() const;					// Second quantization axis of the spin
			arma::vec QuantizationAxis3() const;					// Third quantization axis of the spin
			arma::mat QuantizationAxes() const;						// Transformation matrix consisting of the quantization axes
			
			// Allow access to custom properties to be used for custom tasks
			std::shared_ptr<const MSDParser::ObjectParser> Properties() const;
			
			// Spin matrices
			const arma::sp_cx_mat Sx() const;	// Pauli matrices (Sx, Sy, Sz)
			const arma::sp_cx_mat Sy() const;
			const arma::sp_cx_mat Sz() const;
			const arma::sp_cx_mat Sp() const;	// S+ operator
			const arma::sp_cx_mat Sm() const;	// S- operator
			const arma::sp_cx_mat Tx() const;	// Adjusted magnetic moment operators (Tx, Ty, Tz), which are the Pauli matrices transformed by the tensor
			const arma::sp_cx_mat Ty() const;
			const arma::sp_cx_mat Tz() const;
			
			// Public method for creating ActionTargets
			void GetActionTargets(std::vector<RunSection::NamedActionScalar>&, std::vector<RunSection::NamedActionVector>&, std::string);
	};
	
	// Non-member non-friend functions
	bool HasTrajectory(const Spin& _spin);
	bool IsIsotropic(const Spin& _spin);	// G-tensor isotropy
	
	// Define alias for spin-pointers
	using spin_ptr = std::shared_ptr<Spin>;
}

#endif
