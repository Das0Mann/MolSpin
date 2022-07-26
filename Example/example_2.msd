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
		field = "0.03 0.03 0.025";	      // 0.05 mT along the z-axis
		spins = electron1,electron2;
                tau_c = corr_time;   // Correlation time
                g = 100;                        // Correlation function amplitude
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
		tensor = isotropic(0.4399);	// Isotropic hyperfine tensor of 0.5 mT
		tau_c = corr_time;		                        // Correlation time
		g = 100;			                        // Correlation function amplitude
	}

//	        Interaction HF2
//        {
//                prefactor = 1e-3;
//                type = Hyperfine;
//                group1 = electron2;                             // Spins in group1 interact with spins in group2
//                group2 = nucleus2;
//                tensor = isotropic(1) + anisotropic("0.5 0.5 2");       // Isotropic hyperfine tensor of 0.5 mT
//                tau_c = corr_time;                                   // Correlation time
//                g = 100;                                                // Correlation function amplitude
//        }


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
		steps = 1;	// Make 90 calculation steps (i.e. obtain 90 data points)
	}
}
// -------------------------------------------------------------
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
	
       // BWR-Relaxation method for peturbative treatment of hamiltonians

        Task Method3
        {
                type = sparse-redfield-relaxation;
                logfile = "redfield-relaxation.log";
                datafile = "redfield-relaxation.dat";
		transitionyields = true;

        }	

}
