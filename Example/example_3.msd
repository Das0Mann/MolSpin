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
		spin = 1/2;
	}
	
//	Spin electron2
//	{
//		type = electron;
//	}
	
	// -------------------------
	// Nuclear spins
	// -------------------------
	Spin nucleus1
	{
		spin = 1/2;
		type = nucleus;
	}
	
	// -------------------------
	// Zeeman interaction
	// -------------------------
	Interaction zeeman1
	{
		type = Zeeman;
		field = "0.20 0.30 0.60";	      // 0.05 mT along the z-axis
		tensor = matrix("2.0 2.0 -1.0; 0.0 3.0 0.0; -1.0 2.0 1.0");
		spins = nucleus1;
                tau_c = 0.0098;   // Correlation time
                g = 1;                        // Correlation function amplitud
		coeff = 1;
	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
//	Interaction HF1
//	{
//		type = Hyperfine;
//		group1 = electron1;	                        // Spins in group1 interact with spins in group2
//		group2 = nucleus1;
//		tensor = matrix("40.99 -1.0 -1.0; -1.0 39.99 -1.0; -1.0 -1.0 19.99");	// Isotropic hyperfine tensor of 0.5 mT
//                tau_c = 0.0098;   // Correlation time
//                g = 1;                        // Correlation function amplitud
//		coeff = 1;
//	}

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

	State Identity
	{
	}

        // ---------------------------------------------------------
        // Transitions
        // ---------------------------------------------------------
        Transition Product1
        {
                type = sink;
                source = Singlet; 
                rate = 1.00e-4;
        }

        Transition Product2
        {
                type = sink;
                source = T0;      
                rate = 1.00e-4;
        }

        Transition Product3
        {
                type = sink;
                source = Tm;      
                rate = 1.00e-4;
        }

        Transition Product4
        {
                type = sink;
                source = Tp;      
                rate = 1.00e-4;
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
                steps = 1; //88;
                notifications = details;
        }

}


Run
{
        Task Method3
        {
                type = "redfield-relaxation";
                logfile = "redfield-relaxation.log";
                datafile = "redfield-relaxation.dat";
		transitionyields = true;
        }

}
