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
	// Time-dependent Interactions - Fields
	// -------------------------
	// You can add time-dependent magnetic fields like this.
	// Note that you would need to use a task class that supports time-dependent interactions such as DynamicHS-TimeEvolution.
	
	Interaction linearpolarized       //generate a linearly polarized magnetic field
	{
		type = Zeeman;
		field = "0 0 5e-5";
		spins = electron1, electron2;
		fieldtype = LinearPolarized;
		frequency = 1e-2;
		phase = 0;
        printfield = true; //outputs the file 'linearpolarized.mst'
	}

    // For now, additional time-dependent field interactions are commented out.
    /*
	Interaction circularpolarized      //generate a circularly polarized magnetic field
	{
		type = Zeeman;
		field = "1e-5 0 1e-5";
		spins = electron1, electron2;
		fieldtype = CircularPolarized;
		frequency = 1e-3;
		phase = 1;
		axis = "0 0 1"
		PerpendicularOscillations = false;
	}
    
    Interaction zeeman_BB      //generate a magnetic field modulated by broadband noise
	{
		prefactor = 0.001;
		type = Zeeman;
		field = "0.01 0.0 0.08";
		spins = electron1, electron2;

		fieldtype="broadband";
		minfreq=0.1;
		maxfreq=0.001;
		stdev=0.01;
		components=100;
		randomorientations = false;
		printfield = true;   //outputs the file 'zeeman_BB.mst'
	}
    
    Interaction zeeman_OU      //generate a magnetic field modulated by Ornstein-Uhlenbeck noise
    {
        prefactor = 0.001;
		type = Zeeman;
		field = "0.1 0.0 0.8";
		spins = electron1, electron2;   

        fieldtype = "ougeneral";
        correlationtime = 1e6;
		stdev = 0.3;
		seed = 1;
		timestep = 1; 
    }    
    */

    // -------------------------
	// Time-dependent Interactions - Tensors
	// -------------------------
	// You can add time-dependent interaction tensors like this.

    Interaction HFI_OU     //generate a time-dependent hyperfine interaction modulated by Ornstein-Uhlenbeck noise
	{
	  
		type = DoubleSpin;
        prefactor = 0.001;

		tensor = matrix("-0.099 -0.003 0.000; -0.003 -0.087 0.000; 0.000 0.000 1.757");

		tensortype = "ougeneral";
        correlationtime = 1e6;
		stdev = 1.0;
		seed = 1;
		timestep = 1;

        printtensor = true;    //outputs the file 'HFI_OU.mst'

		group1 = electron1;
		group2 = nucleus1;
	}

    // For now, additional time-dependent tensor interactions are commented out.
    /*
    Interaction HFI_monochromatic  //generate a time-dependent hyperfine interaction modulated by monochromatic noise
    {
        type = DoubleSpin;
        prefactor = 0.001;

		tensor = matrix("-0.099 -0.003 0.000; -0.003 -0.087 0.000; 0.000 0.000 1.757");
		tensortype = "monochromatic";
		frequency = 0.01;
		phase = 0;
		amplitude = 0.1;

        printtensor = false;   

		group1 = electron1;
		group2 = nucleus1;

    }

    Interaction HFI_BB       //generate a time-dependent hyperfine interaction modulated by broadband noise
	{

		type = DoubleSpin;
        prefactor = 0.001;

		tensor = matrix("-0.099 -0.003 0.000; -0.003 -0.087 0.000; 0.000 0.000 1.757");
        tensortype = "broadband";
		minfreq = 0.1;
		maxfreq = 10;
		stdev = 0.5;
		components = 1000;
	
		printtensor = true;   //outputs the file 'HFI_BB.mst'
		seed = 1;
	
		group1 = electron1;	
		group2 = nucleus1;
	}
    */

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

Run
{
	Task Method1
	{
		type = "DynamicHS-Direct-TimeEvo";
		logfile = "log_time_dependent.log";
		datafile = "dat_time_dependent.dat";
		totaltime=8000;
		Timestep=1;
		transitionyields = "false";
		initialstate= "singlet";
	}

}
