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
		prefactor = 0.001;	// Change field units from T to mT/
		type = Zeeman;
		field = "0 0 0.05";	// 0.05 mT along the z-axis
		spins = electron1, electron2;
	}

//		Interaction zeeman_BB
//	{
//		prefactor = 0.001;
//		type = Zeeman;
//		field = "0 0 0.05";
//		spins = electron1, electron2;
//
//		fieldtype="broadband";
//		minfreq=1000;
//		maxfreq=10000;
//		stdev=0.001;
//		components=100;
//	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
//		Interaction HFI_sinat
//	{
//	    prefactor = 0.001;	// Change units from T to mT
//		type = DoubleSpin;
//		tensor = matrix("-0.3, -0.04, -0.1; -0.04,-0.55,0.04; -0.1, 0.04, -0.38");
//		tensortype = sinmat;
//		frequency = 1;
//		phase = 0;
//		group1 = "electron1";	// Spins in group1 interact with spins in group2
//		group2 = "nucleus1";
//		
//	}


	Interaction HF1
	{
	    prefactor = 1;	// Change units from T to mT
		type = DoubleSpin;

		tensor = matrix("-0.0003, -0.00004, -0.0001; -0.00004, -0.00055, 0.00004; -0.0001, 0.00004, -0.00038");
		tensortype = gaussian;
		temperature = 300;
		damping =1e-8;
		restoring = 0.05;
		seed=12;
		timestep=0.01;

		group1 = "electron1";	// Spins in group1 interact with spins in group2
		group2 = "nucleus1";
	}
	
//		Interaction HF2
//	{
//		prefactor = 0.001;	// Change units from T to mT
//		type = Hyperfine;
//		group1 = "electron2";
//		group2 = "nucleus2";
//		tensor = anisotropic(0, 0, 1);	// 1 mT along z-axis
//	}
	
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
		notifications = details;
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
	// Calculate quantum yields using Hilbert-space formalism - general method
	Task Method1
	{
		type = DynamicHS-Direct-TimeEvo;
		logfile = "log_test.log";
		datafile = "dat_test.dat";
		totaltime=10;
		Timestep=0.01;
	}

}
