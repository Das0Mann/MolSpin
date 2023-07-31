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
	}
	
	Spin RPElectron2
	{
		type = electron;
		spin = 1/2;
	}

	Spin FADHBeta1
	{
		spin = 1/2;
		tensor = isotropic("0.4070");
	}

	Spin FADH6
	{
		spin = 1/2;
		tensor = matrix("-0.2569 -0.1273 0.0;-0.1273 -0.4711 0.0;0.0 0.0 -0.4336");
	}

	Spin FADN5
	{
		spin = 1;
		tensor = matrix("-0.0989 0.0039 0.0;0.0039 -0.0881 0.0;0.0 0.0 1.7569");
	}

	Spin FADN10
	{
		spin = 1;
		tensor = matrix("-0.0190 -0.0048 0.0;-0.0048 -0.0196 0.0;0.0 0.0 0.6046");
	}

	// -------------------------
	// Zeeman interaction
	// -------------------------
	
	Interaction zeeman
	{
		prefactor = 0.001;
		type = zeeman;
		field = "0.0 0.0 0.05";
		group1 = RPElectron1,RPElectron2;
	}	

	// -------------------------
	// Hyperfine interactions
	// -------------------------
	Interaction radical1hyperfineFADHBeta1
	{
		prefactor = 0.001;
		type = hyperfine;
		group1 = RPElectron1;
		group2 = FADHBeta1;
		tensor = isotropic(1);
	} 

        Interaction radical1hyperfineFADH6
        {
		prefactor = 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = FADH6;
		tensor = isotropic(1);
        }

        Interaction radical1hyperfineFADN5
        {
		prefactor = 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = FADN5;
		tensor = isotropic(1);
        }
 

        Interaction radical1hyperfineFADN10
        {
		prefactor = 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = FADN10;
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
                rate = 1e-05;
        }

       Transition Product2
        {
                type = sink;
                source = T0;
                rate = 1e-05;
        }

       Transition Product3
        {
                type = sink;
                source = Tp;
                rate = 1e-05;
        }

       Transition Product4
        {
                type = sink;
                source = Tm;
                rate = 1e-05;
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
	// ---------------------------------------------------------
	// General settings
	// ---------------------------------------------------------
	Settings general
	{
		steps = 180;
	}

	// ---------------------------------------------------------
	// Actions
	// ---------------------------------------------------------
	Action scan
	{
		type = rotatevector;
		vector = system1.zeeman.field;
		axis = "0 1 0"; // Rotate the field around the y-axis
		value = 1;
	}

	// ---------------------------------------------------------
	// Outputs objects
	// ---------------------------------------------------------
	Output orientation
	{
		type = vectorangle;
		vector = system1.zeeman.field;
		reference = "0 0 1"; // Get angle relative to z-axis
	}
}

// -------------------------------------------------------------
Run
{
	Task Method1
	{
                type = StaticSS;
                logfile = "example_staticss.log";
                datafile = "example_staticss.dat";
                transitionyields = true;
	}
}

