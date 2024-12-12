// -------------------------------------------------------------
// MolSpin Input file example
// See the user manual for more information
// Find it at www.molspin.eu
// -------------------------------------------------------------
SpinSystem system1
{
	// -------------------------
	// Radical pair electrons
	// -------------------------
	Spin electron1
	{
		type = electron; 		// No need to specify spin then, and g-factor of 2 is default
	}
	
	Spin electron2
	{
		type = electron;
	}
	
	// -------------------------
	// Nuclear spins
	// -------------------------
	Spin nucleus1
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic(1);	// Here we don't want a g-factor of 2
	}
	
	Spin nucleus2
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic(1);
	}
	
	
	// -------------------------
	// Oscillating fields
	// -------------------------
	// You can add oscillating magnetic fields like this.
	// Note that you would need to use a task class that supports time-dependent interactions such as DynamicHS-TimeEvolution.
	// For now, these time-dependent interactions are commented out.
	/*Interaction linearpolarized
	{
		type = Zeeman;
		field = "0 0 5e-5";
		spins = electron1, electron2;
		fieldtype = LinearPolarized;
		frequency = 1e-1;
		phase = 0;
	}
	
	Interaction circularpolarized
	{
		type = Zeeman;
		field = "1e-5 0 1e-5";
		spins = electron1, electron2;
		fieldtype = CircularPolarized;
		frequency = 1e-3;
		phase = 1;
		axis = "0 0 1"
		PerpendicularOscillations = false;
	}*/
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
	Interaction HF1
	{
		prefactor = 0.001;	// Change units from T to mT
		tensortype = SinMat;
		type = Hyperfine;
		frequency = 1e-3;
		phase = 1;
		group1 = "electron1";	// Spins in group1 interact with spins in group2
		group2 = "nucleus1";
	}
	
	Interaction HF2
	{
		prefactor = 0.001;	// Change units from T to mT
		type = Hyperfine;
		group1 = "electron2";
		group2 = "nucleus2";
		tensor = anisotropic(0, 0, 1);	// 1 mT along z-axis
	}
	
	// -------------------------
	// Spin States
	// -------------------------
	// Note that spin states should not be normalized here - that is taken care of by MolSpin.
	// Also note that the spin states of the nuclei are not specified here, which is fine since the state objects
	// are used to produce subspace projection operators, and to set the initial density operator.
	// The initial density operator will then be a mixed ensemble of the nuclear spin states.
	State Singlet
	{
		spins(electron1,electron2) = |1/2,-1/2> - |-1/2,1/2>;
	}
	
	State T0
	{
		spins(electron1,electron2) = |1/2,-1/2> + |-1/2,1/2>;
	}
	
	State Tp
	{
		spin(electron1) = |1/2>;
		spin(electron2) = |1/2>;
	}
	
	State Tm
	{
		spin(electron1) = |-1/2>;
		spin(electron2) = |-1/2>;
	}
	
	State Identity
	{
		// An empty state provides an identity projection
	}
	
	// -------------------------
	// Transitions
	// -------------------------
	// Spin-dependent decays are commented out for now - note that not all task classes can handle spin-dependent decay
	/*Transition singlet_decay
	{
		rate = 1e-3;	// 0.001 /ns, i.e. a lifetime of a microsecond
		source = Singlet;
	}
	
	// Decay from the triplet states must be applied to each triplet state individually
	Transition t0_decay
	{
		rate = 1e-4;
		source = T0;
	}
	
	Transition tp_decay
	{
		rate = 1e-4;
		source = Tp;
	}
	
	Transition tm_decay
	{
		rate = 1e-4;
		source = Tm;
	}*/
	
	// Spin-independent decay can be added using the Identity state defined above
	Transition spinindependent_decay
	{
		rate = 1e-4;
		source = Identity;
	}
	
	// -------------------------
	// Spin system properties
	// -------------------------
	Properties properties
	{
		initialstate = Singlet;
	}
}
// -------------------------------------------------------------
Settings
{
	Settings general
	{
		steps = 1;	// Make 90 calculation steps (i.e. obtain 90 data points)
	}
	
	// Actions are triggered after each calculation step (unless the 'first', 'last' or 'period' parameters are set)
	Action scan_field
	{
		type = RotateVector;
		vector = system1.zeeman1.field;
		axis = "0 1 0";
		value = 1;	// Rotate 1 degree about the axis
	}
	
}

Run
{
	// Calculate quantum yields using Liouville-space formalism - general method
	Task Method1
	{
		type = DynamicHS-TimeEvolution;
		logfile = "log_test.log";
		datafile = "dat_test.dat";
		totaltime=10;
	}

}
