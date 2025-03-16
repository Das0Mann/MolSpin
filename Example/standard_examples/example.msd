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
	// Zeeman interaction
	// -------------------------
	Interaction zeeman1
	{
		prefactor = 0.001;	// Change field units from T to mT
		type = Zeeman;
		field = "0 0 0.05";	// 0.05 mT along the z-axis
		spins = electron1, electron2;
	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
	Interaction HF1
	{
		prefactor = 0.001;	// Change units from T to mT
		type = Hyperfine;
		group1 = "electron1";	// Spins in group1 interact with spins in group2
		group2 = "nucleus1";
		tensor = isotropic(0.5);	// Isotropic hyperfine tensor of 0.5 mT
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
		steps = 90;	// Make 90 calculation steps (i.e. obtain 90 data points)
	}
	
	// Actions are triggered after each calculation step (unless the 'first', 'last' or 'period' parameters are set)
	Action scan_field
	{
		type = RotateVector;
		vector = system1.zeeman1.field;
		axis = "0 1 0";
		value = 1;	// Rotate 1 degree about the axis
	}
	
	// Get the angle between the Zeeman-field and the z-axis
	Output Orientation
	{
		type = vectorangle;
		vector = system1.zeeman1.field;
		reference = "0 0 1";
	}
	
	// Output the components (x, y, z) of the Zeeman-field
	Output BField
	{
		type = components;
		vector = system1.zeeman1.field;
	}
}

Run
{
	// Calculate quantum yields using Liouville-space formalism - general method
	Task Method1
	{
		type = StaticSS;
		logfile = "example_staticss.log";
		datafile = "example_staticss.dat";
	}
	
	// Faster method for calculation of the quantum yields (not Liouville space), but does not work with spin-dependent decays
	Task Method2
	{
		type = StaticHS-SymmetricDecay;
		logfile = "example_statichs_symdec.log";
		datafile = "example_statichs_symdec.dat";
	}
	
	// Fastest quantum yield method, but does not work with spin-dependent decays or coupled electronic spins (i.e. exchange or dipole-dipole interactions)
	Task Method4
	{
		type = RP-SymmetricUncoupled;
		logfile = "example_rpsymuncoupled.log";
		datafile = "example_rpsymuncoupled.dat";
	}
	
	// Calculate eigenvalues and resonance frequencies
	Task Method3
	{
		type = Eigenvalues;
		referencestates = Singlet, T0, Tp, Tm;	// Project eigenstates onto these spin states
		spinlist = electron1, electron2;	// Get transition matrix elements for these spins (resonance analysis)
		//eigenvectors = true;	// Print eigenvectors in the log
		//hamiltonian = true;	// Print full Hamiltonian matrices in the log
		resonances = true;	// Resonance frequency analysis
		separatereal = true;	// Separate real and imaginary parts of output
		//initialtime = 0;
		//totaltime = 100;
		//timestep = 1;
		logfile = "example_eigenvalues.log";	// The log contains resonance frequency data
		datafile = "example_eigenvalues.dat";	// Eigenvalues and eigenstate projections go here
	}
}
