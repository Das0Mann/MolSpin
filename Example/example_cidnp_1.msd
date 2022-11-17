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

	Spin H1
	{
		spin = 1/2;
	}

        Spin C1
        {
                spin = 1/2;
        }

        Spin N1
        {
                spin = 1;
        }


	// -------------------------
	// Zeeman interaction
	// -------------------------
	
	Interaction zeeman1
	{
		type = zeeman;
		field = "0 0 0";
		group1 = RPElectron1;
	}	

        Interaction zeeman2
        {
                type = zeeman;
                field = "0.0 0.0 0";
                spins = RPElectron2;
        }

        Interaction zeemanN1
        {
                type = zeeman;
                field = "0.0 0.0 0";
                spins = N1;
        }


	// -------------------------
	// Hyperfine interactions
	// -------------------------

	Interaction radicalhyperfineH1
	{
		prefactor= 0.001;
		type = hyperfine;
		group1 = RPElectron1;
		group2 = H1;
		tensor = matrix("-0.08 0.00 0.00; 0.00 -0.04 0.00; 0.00 0.00 0.12");
	} 

        Interaction radicalhyperfineC1
        {
                prefactor= 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = C1;
                tensor = matrix("-1.05 0.00 0.00; 0.00 -0.98 0.00; 0.00 0.00 2.03");
        }

        Interaction radicalhyperfineN1
        {
                prefactor= 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = N1;
                tensor = matrix("-0.2 0.00 0.00; 0.00 -0.18 0.00; 0.00 0.00 0.38");
        }

	// -------------------------
        // Dipolar interactions
        // -------------------------

        Interaction radical
        {
                prefactor= 0.001;
                type = doublespin;
                group1 = RPElectron2;
                group2 = RPElectron1;
                tensor = matrix("-0.54 0.0 0.0; 0.0 -0.54 0.0; 0.0 0.0 -0.54");
        }

	// -------------------------
	// Spin States
        // ---------------------------------------------------------
	
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
                rate = 0.01;
        }

        Transition Product2
        {
                type = sink;
                source = Singlet;
                rate = 0.001;
        }


        Transition Product3
        {
                type = sink;
                source = T0;
                rate = 0.001;
        }

        Transition Product4
        {
                type = sink;
                source = Tp;
                rate = 0.001;
        }

        Transition Product5
        {
                type = sink;
                source = Tm;
                rate = 0.001;
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
		steps = 100;
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
	        vector = system1.zeemanN1.field;
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
		nuclei_list=N1;
	}
	
}

