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
	Spin FN5
	{
		spin = 1;
		type = nucleus;
		tensor = isotropic(1);	// Here we don't want a g-factor of 2
	}
	
	Spin WNE1   
	{
		spin = 1;
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
		field = "0.0 0.0 0.05 ";	// 0.05 mT along the z-axis
		spins = electron1, electron2;
	}



	Interaction Dipolar
	{
	   
		type = DoubleSpin;
        prefactor=0.0020023;

		tensor = matrix("0.039 -0.432 0.175; -0.432 -0.279 0.251; 0.175 0.251 0.239");
		tensortype = monochromatic;
		frequency = 0.01;
		phase = 0;
		amplitude = 0.1;
		modulatedistance=true;
	
        printtensor=true;

		group1 = electron1;	// Spins in group1 interact with spins in group2
		group2 = electron2;
	}
    
    Interaction HFI_FN5
	{
	  
		type = DoubleSpin;
        prefactor=0.001;

		tensor = matrix("-0.099 -0.003 0.000; -0.003 -0.087 0.000; 0.000 0.000 1.757");
		tensortype = monochromatic;
		frequency = 0.01;
		phase = 0;
		amplitude = 0.1;
		modulatedistance=false;

        printtensor=false;

		group1 = electron1;
		group2 = FN5;
	}


    Interaction HFI_WNE1
	{
	  
		type = DoubleSpin;
        prefactor=0.001;

		tensor = matrix("-0.053 0.059 -0.046; 0.059  0.564 -0.565; -0.046 -0.565 0.453");
		tensortype = monochromatic;
		frequency = 0.01;
		phase = 0;
		amplitude = 0.1;
		modulatedistance=false;

        printtensor=false;

		group1 = electron1;
		group2 = WNE1;
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
	
	//State Identity
	//{
	//	// An empty state provides an identity projection
	//}
	
	// -------------------------
	// Transitions
	// -------------------------
	// Spin-independent decay can be added using the Identity state defined above
	Transition Product1{type = sink;source = Singlet;rate = 0.001;}
    Transition Product2{type = sink;source = T0;rate = 0.001;}
    Transition Product3{type = sink;source = Tp;rate = 0.001;}
    Transition Product4{type = sink;source = Tm;rate = 0.001;}

	// -------------------------
	// Spin system properties
	// -------------------------
	Properties properties
	{
		initialstate = Singlet;
	}

}

Settings{Settings general {steps = 2; notifications = details;}
        Action scan{
                type = rotatevector;
                vector = system1.zeeman1.field;
                axis = "0 1 0";
                value = 180;
        }
}


Run
{
	// Calculate quantum yields using Hilbert-space formalism - general method
	Task Method1
	{
		type = "DynamicHS-Direct-Yields";
		logfile = "log_tdMono_test.log";
		datafile = "dat_tdMono_test.dat";
		totaltime=5000;
		Timestep=1;
		transitionyields = "false";
		initialstate= "singlet";
	}

}
