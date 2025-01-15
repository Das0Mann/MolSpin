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
	
	
	Interaction zeeman_BB
	{
		prefactor = 0.01;
		type = Zeeman;
		spins = electron1, electron2;

		
        trajectory = "FieldBBNoise.mst";
	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------

	Interaction HF2
	{
		prefactor = 0.001;	// Change units from T to mT
		type = Hyperfine;
		group1 = "electron2";
		group2 = "nucleus2";
		tensor = anisotropic(2, 0, 1);	// 1 mT along z-axis
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
	

	// Spin-independent decay can be added using the Identity state defined above
	Transition spinindependent_decay
	{
		rate = 0;
		source = Singlet;
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

Run
{
	// Calculate quantum yields using Hilbert-space formalism - general method
	Task Method1
	{
		type = "DynamicHS-Direct-TimeEvo";
		logfile = "log_BBfield_mst_test.log";
		datafile = "dat_BBfield_mst_test.dat";
		totaltime=50;
		Timestep=1;
		transitionyields = "false";
		initialstate= "singlet";
	}

}
