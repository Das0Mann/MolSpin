// -------------------------------------------------------------
// MolSpin Input file example
// See the user manual for more information
// Find it at www.molspin.eu
// -------------------------------------------------------------
SpinSystem system1
{
	// ---------------------------------------------------------
	// Spins
	// ---------------------------------------------------------

	Spin RPElectron1
	{
		type = electron;
		spin = 1/2;
		tensor = isotropic(2.0035);
	}
	
	Spin RPElectron2
	{
		type = electron;
		spin = 1/2;
		tensor = isotropic(2.002);
	}

	Spin C1
	{
		type = nucleus;
		spin = 1/2;
	}

	// -------------------------
	// Zeeman interaction
	// -------------------------
	
	Interaction zeeman1
	{
		type = zeeman;
		field = "0 0 1";
		group1 = RPElectron1;
	}	

        Interaction zeeman2
        {
                type = zeeman;
                field = "0.0 0.0 1";
                spins = RPElectron2;
        }

        Interaction zeemanC1
        {
                type = zeeman;
                field = "0.0 0.0 1";
                spins = C1;
		CommonPrefactor = true;
		Prefactor = -0.000382102;
        }

	// -------------------------
	// Hyperfine interactions
	// -------------------------

	Interaction radicalhyperfineC1
	{
		type = hyperfine;
		group1 = RPElectron1;
		group2 = C1;
		tensor = matrix("-0.099 -0.003 0.000;-0.003 -0.087 0.000;0.000 0.000 1.757");
		commonprefactor = true;
		prefactor= 0.001;
	} 

	// -------------------------
        // Dipolar interactions
        // -------------------------

        Interaction dipolar
        {
                type = doublespin;
                group1 = RPElectron1;
                group2 = RPElectron2;
		tensor=matrix("-0.002 0.0 0.0;0.0 -0.002 0.0;0.0 0.0 -0.002");
		commonprefactor = true;	
                prefactor= 0.001;
        }

	// -------------------------
	// Spin States
        // -------------------------
	
        State Singlet	// |S>
	{
		spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;
	}
	
	State T0	// |T0>
	{
		spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;
	}
	
	State Tp	// |T+>
	{
		spin(RPElectron2) = |1/2>;
		spin(RPElectron1) = |1/2>;
	}
	
	State Tm	// |T->
	{
		spin(RPElectron2) = |-1/2>;
		spin(RPElectron1) = |-1/2>;
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
                rate = 0.015;
        }

        Transition Product2
        {
                type = sink;
                source = T0;
                rate = 0.015;
        }

        Transition Product3
        {
                type = sink;
                source = Tp;
                rate = 0.015;
        }

        Transition Product4
        {
                type = sink;
                source = Tm;
                rate = 0.015;
        }

        // -------------------------
	// Spin system properties
	// -------------------------

	Properties properties
	{
		initialstate = T0;
	}
}
// -------------------------------------------------------------
Settings
{
	Settings general
	{
		steps = 200;
		notifications = details;
	}

        Action field_strength1
        {
	        type = AddVector;
	        vector = system1.zeeman1.field;
	        direction = "0 0 1";
        	value = 0.1;
        }

        Action field_strength2
        {
	        type = AddVector;
	        vector = system1.zeeman2.field;
	        direction = "0 0 1";
        	value = 0.1;
        }

        Action field_strength2
        {
	        type = AddVector;
	        vector = system1.zeemanC1.field;
	        direction = "0 0 1";
        	value = 0.1;
        }
}
// -------------------------------------------------------------
Run
{
	Task Method1
	{
                type = StaticSS-CIDNP;
                logfile = "example_staticsscidnp.log";
                datafile = "example_staticsscidnp.dat";
                transitionyields = true;
		nuclei_list=C1;
	}	
}

