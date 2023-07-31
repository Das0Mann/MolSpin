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
		tensor = isotropic(2.0023);
	}
	
	Spin RPElectron2
	{
		type = electron;
		spin = 1/2;
		tensor = isotropic(2.0023);
	}

	Spin H1
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
                field = "0 0 50";
                spins = RPElectron1;
        }

	Interaction zeeman2
	{
		type = zeeman;
		field = "0 0 50";
		spins = RPElectron2;
	}
	
	// -------------------------
	// Hyperfine interactions
	// -------------------------
	Interaction radical1hyperfine
	{
		type = hyperfine;
		group1 = RPElectron1;
		group2 = H1;
                tensor = matrix("0.5 0.0 0.0;0.0 0.5 0.0;0.0 0.0 2.0");
		tau_c = 0.0001;   
		g = 1;      
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

        Transition Product4
        {
                type = sink;
                source = T0;
                rate =0.001;
        }

        Transition Product5
        {
                type = sink;
                source = Tp;
                rate =0.001;
        }

        Transition Product6
        {
                type = sink;
                source = Tm;
                rate =0.001;
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
// -------------------------------------------------------------
Run
{
	Task Method2
	{
		type = redfield-relaxation;
                logfile = "redfield-relaxation.log";
                datafile = "redfield-relaxation.dat";
                transitionyields = true;
	
	}

}

