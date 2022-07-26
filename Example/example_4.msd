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
		spin = 5/2;
		type = nucleus;
		tensor = isotropic(1);	// Here we don't want a g-factor of 2
	}
	
	// -------------------------
	// Zeeman interaction
	// -------------------------
	Interaction zeeman1
	{
		prefactor = 1e-3;
		type = Zeeman;
		field = "0.00 0.00 50";	      // 0.05 mT along the z-axis
		spins = electron1;
                tau_c = 0.0098,0.0098,0.0098,0.0098,0.0098,0.0098;   // Correlation time
                g = 0.001,0.001,0.001,0.001,0.001,0.001;                        // Correlation function amplitud
		terms = 0;
		define_ampl = 1;
	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
	Interaction HF1
	{
		prefactor = 1e-3;
		type = Hyperfine;
		group1 = electron1,electron2;	                        // Spins in group1 interact with spins in group2
		group2 = nucleus1;
		tensor = isotropic(43.99);	// Isotropic hyperfine tensor of 0.5 mT
                tau_c = 0.0098,0.0098,0.0098,0.0098,0.0098,0.0098;   // Correlation time
                g = 0.001,0.01,0.1,1,1,1;                        // Correlation function amplitud
		ops = 0;
		define_ampl = 1;
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
	
	Action rotatevector
	{
		type = RotateVector;
		value = 0.01;
		axis = "0 1 0";
		vector = system1.zeeman1.field;
	}

}
// -------------------------------------------------------------
Run
{
	// Calculate quantum yields using Liouville-space formalism - general method
//	Task Method1
//	{
//		type = "staticss-timeevolution";
//		logfile = "example_timeevolution.log";
//		datafile = "example_timeevolution.dat";
//		Timestep = 0.1;
//		TotalTime = 250;
//	}
	
       // BWR-Relaxation method for peturbative treatment of hamiltonians

        Task Method2
        {
                type = "redfield-relaxation-timeevolution";
                logfile = "redfield-relaxation-timeevolution.log";
                datafile = "redfield-relaxation-timeevolution.dat";
	        Timestep = 0.1;
                TotalTime = 250;
        }	


        Task Method3
        {
                type = "redfield-relaxation";
                logfile = "redfield-relaxation.log";
                datafile = "redfield-relaxation.dat";
		transitionyields = true;
        }

       Task Method4
        {
                type = "sparse-redfield-relaxation";
                logfile = "sparse-redfield-relaxation.log";
                datafile = "sparse-redfield-relaxation.dat";
                transitionyields = true;
        }

        Task Method5
        {
                type = "staticss";
                logfile = "staticSS.log";
                datafile = "staticSS.dat";
                transitionyields = true;
        }

}
