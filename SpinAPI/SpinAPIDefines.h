/////////////////////////////////////////////////////////////////////////
// Defines (SpinAPI Module)
// ------------------
// Definitions used by the SpinAPI Module.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_SpinAPI_Defines
#define MOD_SpinAPI_Defines

namespace SpinAPI
{
	// Used by the Spin class
	enum class SpinType
	{
		Electron,
		Nucleus,
		NotSpecified,
	};

	// Used by the Trajectory class
	enum class InterpolationType
	{
		Stepwise, // Stepwise constant function
		Linear,	  // Linear interpolation
	};

	// Used by the Interaction class
	enum class InteractionType
	{
		Undefined,
		SingleSpin,
		DoubleSpin,
		Exchange,
		Zfs,
	};

	// Used by the Interaction class to determine the time-dependence of the field for SingleSpin interactions
	enum class InteractionFieldType
	{
		Static,
		LinearPolarization,	  // Monochromatic linearly polarized radiation, parameters: "frequency", "phase"
		CircularPolarization, // Monochromatic circularly polarized radiation, parameters: "frequency", "phase", "axis"
		BroadbandNoise,
		Trajectory,			  // Time-dependence is given by a trajectory
	};

	enum class InteractionTensorType //PB added
	{
		Static,
		SinMat, //example 
		GaussianNoise,
		BroadbandNoise,
		Trajectory,			   // Time-dependence is given by a trajectory
	};

	// Used by the Transition class
	enum class TransitionType
	{
		Source,
		Sink,
	};

	// The types of supported reaction operators
	enum class ReactionOperatorType
	{
		Unspecified, // Used by the transition class - SpinAPI::SpinSpace will use the reaction operator type assigned to it if a transition has an unspecified reaction operator type
		Haberkorn,
		Lindblad,
	};

	// The types of special operators defined in SpinAPI::Operator objects
	enum class OperatorType
	{
		Unspecified,
		RelaxationLindblad, // Single-spin operator, i.e. uses Sx, Sy and Sz operators of the specified spins
		RelaxationDephasing,
		RelaxationRandomFields,
		RelaxationT1,
		RelaxationT2,
	};

	// The types of special operators defined in SpinAPI::Operator objects
	enum class PulseType
	{
		Unspecified,
		InstantPulse, // Single-spin operator, i.e. uses Sx, Sy and Sz operators of the specified spins
		LongPulse,
		LongPulseStaticField,
		ShapedPulse,
	};

	// Types of standard outputs based on ActionTargets, to be used to by the StandardOutput class
	enum class StandardOutputType
	{
		VectorXYZ,
		VectorAngle,
		VectorLength,
		VectorDot,
		Scalar,
		Undefined,
	};
}

#endif
