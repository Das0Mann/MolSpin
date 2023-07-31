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
		tensor = matrix("0.00024374807218487007 -3.080644368944783e-05 2.9091859084380416e-05;-3.080644368944782e-05 0.00013035721866991337 -5.754234641385961e-06;2.909185908438041e-05 -5.75423464138596e-06 0.00011702083253180584");
	}

	Spin FADH6
	{
		spin = 1/2;
		tensor = matrix("-0.0001833333306189016 3.8247743672635986e-05 -1.5375003551147258e-06;3.8247743672635986e-05 -0.0004576429525705895 -2.353081542319135e-07;-1.537500355114725e-06 -2.3530815423191244e-07 -0.0004165071949152782");
	}

	Spin FADN5
	{
		spin = 1;
		tensor = matrix("-1.09302886e-05 -2.81218002e-06 -6.96255148e-05;-2.81218002e-06 -2.22239306e-05 1.15009400e-05;-6.96255148e-05 1.15009400e-05 0.00179211");
	}

	Spin FADN10
	{
		spin = 1;
		tensor = matrix("1.72116133e-05 -1.66650515e-06 -2.60628815e-06;-1.66650515e-06 6.5802297e-06 -5.54980951e-06;-2.60628815e-06 -5.54980951e-06 0.00063954");
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
//		prefactor = 0.001;
		type = hyperfine;
		group1 = RPElectron1;
		group2 = FADHBeta1;
	} 

        Interaction radical1hyperfineFADH6
        {
//		prefactor = 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = FADH6;
        }

        Interaction radical1hyperfineFADN5
        {
//		prefactor = 0.001;
                type = hyperfine;
                group1 = RPElectron1;
                group2 = FADN5;
        }

        Interaction radical1hyperfineFADN10
        {
//		prefactor = 0.001;
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
                rate = 0.001;
        }

       Transition Product2
        {
                type = sink;
                source = T0;
                rate = 0.001;
        }

       Transition Product3
        {
                type = sink;
                source = Tp;
                rate = 0.001;
        }

       Transition Product4
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

