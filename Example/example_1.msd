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
	
	// -------------------------
	// Zeeman interaction
	// -------------------------
	Interaction zeeman1
	{
		prefactor = 0.001;	      // Change field units from T to mT
		type = Zeeman;
		field = "0 0 0";	      // 0.05 mT along the z-axis
		spins = electron1,electron2;
                tau_c = 0.001;   // Correlation time
                g = 1;                        // Correlation function amplitude
	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
	Interaction HF1
	{
		prefactor = 1e-3;	
		type = Hyperfine;
		group1 = electron1;	                        // Spins in group1 interact with spins in group2
		group2 = nucleus1;
		tensor = isotropic(1) + anisotropic("0.5 0.5 2");	// Isotropic hyperfine tensor of 0.5 mT
		tau_c = 0.001;		                        // Correlation time
		g = 1;			                        // Correlation function amplitude
	}

	        Interaction HF2
        {
                prefactor = 1e-3;
                type = Hyperfine;
                group1 = electron1;                             // Spins in group1 interact with spins in group2
                group2 = nucleus2;
                tensor = isotropic(1) + anisotropic("0.5 0.5 2");       // Isotropic hyperfine tensor of 0.5 mT
                tau_c = 0.001;                                   // Correlation time
                g = 1;                                                // Correlation function amplitude
        }


	// -------------------------
	// Spin States
	// -------------------------
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

	State mp
	{
                spin(electron1) = |1/2>;
                spin(electron2) = |-1/2>;
	}
	
	State pm
	{
                spin(electron1) = |1/2>;
                spin(electron2) = |1/2>;
	}
	
	State Identity
	{
	}

        // ---------------------------------------------------------
        // Transitions
        // ---------------------------------------------------------
        Transition Product1
        {
                type = sink;
                source = Identity;      // spin-independent reaction
                rate = 1.00e+6;
        }

	// -------------------------
	// Spin system properties
	// -------------------------
	Properties properties
	{
		initialstate = Tm;
	}
}
// -------------------------------------------------------------
Settings
{
	Settings general
	{
		steps = 1;	// Make 90 calculation steps (i.e. obtain 90 data points)
	}
	
	Action scan
	{
		type = rotatevector;
		vector = system1.zeeman1.field;
		axis = "0 1 0";
		value= 0.5;
	}

	Output orientation
	{
		type = vectorangle;
		vector = system1.zeeman1.field;
		reference = "0 0 1"
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
		transitionyields = true;
	}
	
	Task main
	{ 
		type = RP-SymmetricUncoupled; // Specify task class
		logfile = "example_RP-SymmetricUncoupled.log";
		datafile = "example_RP-SymmetricUncoupled.dat"; 
	}

	
       // BWR-Relaxation method for peturbative treatment of hamiltonians

        Task Method3
        {
                type = redfield-relaxation;
                logfile = "redfield-relaxation.log";
                datafile = "redfield-relaxation.dat";
		transitionyields = true;

        }	

}
