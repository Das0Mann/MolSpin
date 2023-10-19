/////////////////////////////////////////////////////////////////////////
// SpinSpace class (SpinAPI Module)
// ------------------
// Helper class that defines operators/matrices on the space spanned by
// a collection of spins. E.g. the Sz operator in the total Hilbert space.
// 
// The SpinSpace class contains methods to produce operators in various
// formats (dense, sparse, Liouville-space-dense, Liouville-space-sparse),
// but for distributed matrix formats use the DistributedSpinSpace class	// TODO: Create DistributedSpinSpace class
// instead (can be constructed from a SpinSpace).
// 
// This is meant as a helper class that can be used in the calculation
// tasks, but it is not necessary to use it - e.g. semi-classical methods
// will not be using it. The SpinSpace class is not a friend class of
// any other classes, so every task can be implemented without the use of
// SpinSpace - the SpinSpace class is just here to make life easier.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_SpinSpace
#define MOD_SpinAPI_SpinSpace

#include <vector>
#include <memory>
#include <armadillo>
#include "SpinAPIDefines.h"
#include "SpinAPIfwd.h"

namespace SpinAPI
{
	class SpinSpace
	{
		private:
			// Implementation
			bool useSuperspace;
			std::vector<spin_ptr> spins;
			std::vector<interaction_ptr> interactions;
			std::vector<transition_ptr> transitions;
			double time;					// The current time to use
			unsigned int trajectoryStep;	// The current trajectory step to use
			bool useTrajectoryStep;			// Set to true if trajectories should be used instead of time, where available
			ReactionOperatorType reactionOperators;
			
		public:
			// Constructors / Destructors
			SpinSpace();										// Normal constructors
			explicit SpinSpace(const spin_ptr&);
			explicit SpinSpace(const std::vector<spin_ptr>&);
			explicit SpinSpace(const std::shared_ptr<SpinSystem>&);
			explicit SpinSpace(const SpinSystem&);
			SpinSpace(const SpinSpace&);						// Copy-constructor
			~SpinSpace();										// Destructor
			
			// Operators
			const SpinSpace& operator=(const SpinSpace&);		// Copy-assignment
			
			// ------------------------------------------------
			// Spin management (SpinSpace_management.cpp)
			// ------------------------------------------------
			// Adding spins to the space
			bool Add(const spin_ptr&);							// Add the spin to the space
			bool Add(const std::vector<spin_ptr>&);				// Add all spins in the list to space (spins will not be duplicated)
			bool Add(const std::shared_ptr<SpinSystem>&);		// Add all spins and interactions from the SpinSystem (spins will not be duplicated)
			bool Add(const SpinSystem&);
			bool CompleteSpace(const state_ptr&);				// Add missing spins (if any) to make the space comeplete with respect to the state
			bool CompleteSpace(const interaction_ptr&);			// Add missing spins (if any) to make the space comeplete with respect to the interaction
			
			// Removing spins from the space
			bool Remove(const spin_ptr&);						// Removes the spin from the space
			bool Remove(const std::vector<spin_ptr>&);			// Removes all spins in the list from the space
			bool Remove(const std::shared_ptr<SpinSystem>&);	// Removes all spins in the SpinSystem from the space
			bool Remove(const SpinSystem&);
			bool RemoveSubspace(const spin_ptr&, const state_ptr&);					// Removes the spin and all spins coupled/entangled with it in the given state from the space
			bool RemoveSubspace(const spin_ptr&, const interaction_ptr&);			// Removes the spin and all spins interacting with it in the given interaction from the space
			void ClearSpins();														// Removes all spins from the spin space
			
			// Checking whether the space contains a range of spins
			bool Contains(const spin_ptr&) const;									// Checks single spin
			bool Contains(const std::vector<spin_ptr>&) const;						// Checks whether all spins in range are contained in the space
			bool Contains(const std::shared_ptr<SpinSystem>&) const;				// Checks whether all spins in the SpinSystem are contained in the space
			bool Contains(const SpinSystem&) const;
			bool ContainsSubspace(const spin_ptr&, const state_ptr&) const;			// Checks whether all spins coupled/entangled with the spin are contained in the space
			bool ContainsSubspace(const spin_ptr&, const interaction_ptr&) const;	// Checks whether all spins interacting with the spin are contained in the space
			
			// ------------------------------------------------
			// Interaction management (SpinSpace_management.cpp)
			// ------------------------------------------------
			// Adding interactions to the space
			bool Add(const interaction_ptr&);							// Add the interaction to the space
			bool Add(const std::vector<interaction_ptr>&);				// Add all interactions in the list to space (interactions will not be duplicated)
			bool Remove(const interaction_ptr&);						// Removes the interaction from the space
			bool Remove(const std::vector<interaction_ptr>&);			// Removes all interactions in the list from the space
			bool Contains(const interaction_ptr&) const;				// Checks single interaction
			bool Contains(const std::vector<interaction_ptr>&) const;	// Checks whether all interactions in range are contained in the space
			void ClearInteractions();									// Removes all interactions from the spin space
			
			// ------------------------------------------------
			// Transition management (SpinSpace_management.cpp)
			// ------------------------------------------------
			// Adding transitions to the space
			bool Add(const transition_ptr&);							// Add the transition to the space
			bool Add(const std::vector<transition_ptr>&);				// Add all transitions in the list to space (transitions will not be duplicated)
			bool Remove(const transition_ptr&);							// Removes the transition from the space
			bool Remove(const std::vector<transition_ptr>&);			// Removes all transitions in the list from the space
			bool Contains(const transition_ptr&) const;					// Checks single transition
			bool Contains(const std::vector<transition_ptr>&) const;	// Checks whether all transitions in range are contained in the space
			void ClearTransitions();									// Removes all transitions from the spin space
			
			// ------------------------------------------------
			// Spin state representations in the Hilbert space (SpinSpace_states.cpp)
			// ------------------------------------------------
			arma::cx_mat GetSingleSpinState(const spin_ptr&, int) const;	// Returns the projection matrix
			bool GetState(const CompleteState&, arma::cx_vec&, bool _useFullBasis = true) const;			// Vector representing the state
			bool GetState(const state_ptr&, arma::cx_vec&) const;			// Vector representing the state
			bool GetState(const state_ptr&, arma::cx_mat&) const;			// Projection operator onto the state (dense matrix)
			bool GetState(const state_ptr&, arma::sp_cx_mat&) const;		// Projection operator onto the state (sparse matrix)
			
			// ------------------------------------------------
			// Operators in the spin space (SpinSpace_operators.cpp)
			// ------------------------------------------------
			bool CreateOperator(const arma::cx_mat&, const spin_ptr&, arma::cx_mat&) const;					// Creates an operator in the Hilbert space. The matrix must be square with the multiplicity of the spin as its dimension
			bool OperatorToSuperspace(const arma::cx_mat&, arma::cx_vec&) const;							// Converts the operator to a vector in the superspace
			bool OperatorFromSuperspace(const arma::cx_vec&, arma::cx_mat&) const;							// Converts a superspace vector back to a Hilbert space operator
			bool SuperoperatorFromOperators(const arma::cx_mat&, const arma::cx_mat&, arma::cx_mat&) const;	// Creates a superspace operator from the two given (left-side and right-side) operators
			bool SuperoperatorFromLeftOperator(const arma::cx_mat&, arma::cx_mat&) const;					// Assumes right-side operator is identity (more efficient than passing the identity to SuperoperatorFromOperators)
			bool SuperoperatorFromRightOperator(const arma::cx_mat&, arma::cx_mat&) const;					// Assumes left-side operator is identity (more efficient than passing the identity to SuperoperatorFromOperators)

			// Sparse versions
			bool CreateOperator(const arma::sp_cx_mat&, const spin_ptr&, arma::sp_cx_mat&) const;
			bool OperatorToSuperspace(const arma::sp_cx_mat&, arma::cx_vec&) const;
			bool OperatorFromSuperspace(const arma::cx_vec&, arma::sp_cx_mat&) const;
			bool SuperoperatorFromOperators(const arma::sp_cx_mat&, const arma::sp_cx_mat&, arma::sp_cx_mat&) const;
			bool SuperoperatorFromLeftOperator(const arma::sp_cx_mat&, arma::sp_cx_mat&) const;
			bool SuperoperatorFromRightOperator(const arma::sp_cx_mat&, arma::sp_cx_mat&) const;
			
			// Re-ordering of the spins, used by the GetState methods when working with entangled states
			// TODO: Consider making non-member non-friend functions, or making such equivalents
			bool ReorderBasis(arma::cx_vec&, const std::vector<spin_ptr>&) const;
			bool ReorderBasis(arma::cx_mat&, const std::vector<spin_ptr>&) const;
			bool ReorderBasis(arma::sp_cx_mat&, const std::vector<spin_ptr>&) const;
			bool ReorderBasis(arma::cx_vec&, const std::vector<spin_ptr>&, const std::vector<spin_ptr>&) const;
			bool ReorderBasis(arma::cx_mat&, const std::vector<spin_ptr>&, const std::vector<spin_ptr>&) const;
			bool ReorderBasis(arma::sp_cx_mat&, const std::vector<spin_ptr>&, const std::vector<spin_ptr>&) const;
			bool ReorderingOperator(arma::sp_cx_mat&, const std::vector<spin_ptr>&, const std::vector<spin_ptr>&) const;	// Returns the reordering operator itself, used by the previous methods
			
			// ------------------------------------------------
			// Spherical tensors (SpinSpace_operators.cpp)
			// ------------------------------------------------

			// Rank 1			
			bool Rk1SphericalTensorT0(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const;
			bool Rk1SphericalTensorTp1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const;
			bool Rk1SphericalTensorTm1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const;
			
			//Linear interactions rank 0 & 2
			bool LRk0TensorT0(const spin_ptr& _spin1,const arma::cx_vec& _field, arma::cx_mat& _out) const;
			bool LRk0TensorT0(const spin_ptr& _spin1,const arma::cx_vec& _field, arma::sp_cx_mat& _out) const;
			bool LRk2SphericalTensorT0(const spin_ptr&, const arma::cx_vec& _field, arma::cx_mat&) const;		// T(m=0)
			bool LRk2SphericalTensorT0(const spin_ptr&, const arma::cx_vec& _field, arma::sp_cx_mat&) const;	// T(m=0), sparse
			bool LRk2SphericalTensorTp1(const spin_ptr&, const arma::cx_vec& _field, arma::cx_mat&) const;		// T(m=+1)
			bool LRk2SphericalTensorTp1(const spin_ptr&, const arma::cx_vec& _field, arma::sp_cx_mat&) const;	// T(m=+1), sparse
			bool LRk2SphericalTensorTm1(const spin_ptr&, const arma::cx_vec& _field, arma::cx_mat&) const;		// T(m=-1)
			bool LRk2SphericalTensorTm1(const spin_ptr&, const arma::cx_vec& _field, arma::sp_cx_mat&) const;	// T(m=-1), sparse
			bool LRk2SphericalTensorTp2(const spin_ptr&, const arma::cx_vec& _field, arma::cx_mat&) const;		// T(m=+2)
			bool LRk2SphericalTensorTp2(const spin_ptr&, const arma::cx_vec& _field, arma::sp_cx_mat&) const;	// T(m=+2), sparse
			bool LRk2SphericalTensorTm2(const spin_ptr&, const arma::cx_vec& _field, arma::cx_mat&) const;		// T(m=-2)
			bool LRk2SphericalTensorTm2(const spin_ptr&, const arma::cx_vec& _field, arma::sp_cx_mat&) const;	// T(m=-2), sparse

			//Bilinear interactions rank 0 & 2
			bool BlRk0TensorT0(const spin_ptr& _spin1,const spin_ptr& _spin2 , arma::cx_mat& _out) const;
			bool BlRk0TensorT0(const spin_ptr& _spin1,const spin_ptr& _spin2 , arma::sp_cx_mat& _out) const;
			bool BlRk2SphericalTensorT0(const spin_ptr&, const spin_ptr&, arma::cx_mat&) const;		// T(m=0)
			bool BlRk2SphericalTensorT0(const spin_ptr&, const spin_ptr&, arma::sp_cx_mat&) const;	// T(m=0), sparse
			bool BlRk2SphericalTensorTp1(const spin_ptr&, const spin_ptr&, arma::cx_mat&) const;		// T(m=+1)
			bool BlRk2SphericalTensorTp1(const spin_ptr&, const spin_ptr&, arma::sp_cx_mat&) const;	// T(m=+1), sparse
			bool BlRk2SphericalTensorTm1(const spin_ptr&, const spin_ptr&, arma::cx_mat&) const;		// T(m=-1)
			bool BlRk2SphericalTensorTm1(const spin_ptr&, const spin_ptr&, arma::sp_cx_mat&) const;	// T(m=-1), sparse
			bool BlRk2SphericalTensorTp2(const spin_ptr&, const spin_ptr&, arma::cx_mat&) const;		// T(m=+2)
			bool BlRk2SphericalTensorTp2(const spin_ptr&, const spin_ptr&, arma::sp_cx_mat&) const;	// T(m=+2), sparse
			bool BlRk2SphericalTensorTm2(const spin_ptr&, const spin_ptr&, arma::cx_mat&) const;		// T(m=-2)
			bool BlRk2SphericalTensorTm2(const spin_ptr&, const spin_ptr&, arma::sp_cx_mat&) const;	// T(m=-2), sparse

			// ------------------------------------------------
			// New added functions for wavefucntion formalism and SSE
			// ------------------------------------------------
			
			arma::cx_colvec SUZstate(const int &spinmult, std::mt19937& generator); // returns stochastically determined SU(Z) state  	
			arma::cx_colvec CoherentState(std::vector<SpinAPI::system_ptr>::const_iterator i, std::mt19937& generator);
			arma::cx_mat HighamProp(arma::sp_cx_mat &H, arma::cx_mat& B, const std::complex<double> t, const std::string precision, arma::mat &M);
			arma::mat SelectTaylorDegree(const arma::sp_cx_mat &H, const std::string precision, const int lengthB);
			double normAmEst(const arma::sp_cx_mat &H, double m, std::mt19937& generator);
			arma::cx_colvec KrylovExpmGeneral(const arma::sp_cx_mat &H, const arma::cx_colvec &b, const arma::cx_double dt, int KryDim, int HilbSize);
			arma::cx_colvec KrylovExpmSymm(const arma::sp_cx_mat &H, const arma::cx_colvec &b, const arma::cx_double dt, int KryDim, int HilbSize);						
			void ArnoldiProcess(const arma::sp_cx_mat& H, const arma::cx_colvec& b, arma::cx_mat& KryBasis, arma::cx_mat& Hessen, int KryDim, double& h_mplusone_m);	
			void LanczosProcess(const arma::sp_cx_mat& H, const arma::cx_colvec& b, arma::cx_mat& KryBasis, arma::cx_mat& Hessen, int KryDim, double& h_mplusone_m);

			// ------------------------------------------------
			// Hamiltonian representations in the space (SpinSpace_hamiltonians.cpp)
			// ------------------------------------------------
			bool InteractionOperator(const interaction_ptr&, arma::cx_mat&) const;		// Returns the matrix representation of the interaction on the spin space (dense matrix)
			bool InteractionOperator(const interaction_ptr&, arma::sp_cx_mat&) const;	// Returns the matrix representation of the interaction on the spin space (sparse matrix)
			bool Hamiltonian(arma::cx_mat&) const;				// Total Hamiltonian operator (dense matrix)
			bool Hamiltonian(arma::sp_cx_mat&) const;			// Total Hamiltonian operator (sparse matrix)
			bool StaticHamiltonian(arma::cx_mat&) const;		// Time-independent part of the Hamiltonian operator (dense matrix)
			bool StaticHamiltonian(arma::sp_cx_mat&) const;		// Time-independent part of the Hamiltonian operator (sparse matrix)
			bool DynamicHamiltonian(arma::cx_mat&) const;		// Time-dependent part of the Hamiltonian operator (dense matrix)
			bool DynamicHamiltonian(arma::sp_cx_mat&) const;	// Time-dependent part of the Hamiltonian operator (sparse matrix)
			
			// ------------------------------------------------
			// Transitions/decay operators (SpinSpace_transitions.cpp)
			// ------------------------------------------------
			// NOTE: Only works with Haberkorn operators! TODO: Consider implementations for other reaction operator types.
			// Provides the operator "k/2 * P" in Hilbert space, where P is the projection onto a state
			// Provides the anti-commutator operator "k/2 * {P,.}" in superoperator space
			bool ReactionOperator(const transition_ptr&, arma::cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;		// The last parameter can be used to force the use of a specific reaction operator type,
			bool ReactionOperator(const transition_ptr&, arma::sp_cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;	// but by default the reaction operator type of the transition or the spinspace will be used.
			bool TotalReactionOperator(arma::cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;			// Total reaction operator (dense matrix)
			bool TotalReactionOperator(arma::sp_cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;		// Total reaction operator (sparse matrix)
			bool StaticTotalReactionOperator(arma::cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;		// Time-independent part of the total reaction operator (dense matrix)
			bool StaticTotalReactionOperator(arma::sp_cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;	// Time-independent part of the total reaction operator (sparse matrix)
			bool DynamicTotalReactionOperator(arma::cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;	// Time-dependent part of the total reaction operator (dense matrix)
			bool DynamicTotalReactionOperator(arma::sp_cx_mat&, const ReactionOperatorType& _forcedReactionOperatorType = ReactionOperatorType::Unspecified) const;	// Time-dependent part of the total reaction operator (sparse matrix)
			ReactionOperatorType GetReactionOperatorType() const;		// Returns the reaction operator type used by the SpinSpace (superspace only)
			
			// Methods to create reaction operators in the target spin system (i.e. for creation), where the 'double' describes the amount of source state in the source system
			bool ReactionTargetOperator(const transition_ptr&, double, arma::cx_mat&) const;
			bool ReactionTargetOperator(const transition_ptr&, double, arma::sp_cx_mat&) const;
			
			// ------------------------------------------------
			// Relaxation operators (SpinSpace_relaxation.cpp)
			// ------------------------------------------------
			// NOTE: These operators can only be created in superoperator space!
			bool RelaxationOperator(const operator_ptr&, arma::cx_mat&) const;
			bool RelaxationOperator(const operator_ptr&, arma::sp_cx_mat&) const;
			
			// ------------------------------------------------
			// Other public methods
			// ------------------------------------------------
			unsigned int SpaceDimensions() const;						// Returns the size of the spin space, depending on whether superspace is used or not
			unsigned int HilbertSpaceDimensions() const;				// Returns the size of the spin space Hilbert space
			unsigned int SuperSpaceDimensions() const;					// Returns the size of the spin space super-space
			bool HasTimedependentInteractions() const;
			bool HasTimedependentTransitions() const;
			
			// ------------------------------------------------
			// Settings for the spin space
			// ------------------------------------------------
			bool UseSuperoperatorSpace(bool);							// Determines whether all returned operators, etc. will be in superoperator-/Liouville-space
			bool SetReactionOperatorType(const ReactionOperatorType&);	// Sets the type of reaction operator to be produced - NOTE: Only works in superspace
			bool SetTime(double);										// Set the current time, used to set states from trajectories (provided the trajectories have "time" columns, otherwise first step is used)
			bool SetTrajectoryStep(unsigned int);						// Set the current step to be used in all trajectories (trajectories with too few steps will use last step)
	};
	
	// Non-member non-friend functions
	// ------------------------------------------------
	// Creation operators (SpinSpace_transitions.cpp)
	// ------------------------------------------------
	// The creation operators have the form C = |a><b| where |b> is a state in the source system,
	// and |a> is a state in the target system. Thus this operator transforms state |b> of the source
	// system into state |a> of the target system.
	// NOTE: Creation operators do not contain the rate constant.
	bool CreationOperator(const transition_ptr&, const SpinSpace&, const SpinSpace&, arma::cx_mat&, bool _useSuperoperatorSpace = false);
	bool CreationOperator(const transition_ptr&, const SpinSpace&, const SpinSpace&, arma::sp_cx_mat&, bool _useSuperoperatorSpace = false);
}

#endif
