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

    Spin FN5
    {
        tensor = matrix("-9.949848662047867e-05 -2.8710659792386205e-06 0.0;-2.8710659792386357e-06 -8.749057039785496e-05 0.0;0.0 0.0 0.0017568790352424672");
        spin = 1;
    }

    Spin WCNE1
    {
        tensor = matrix("-5.2960934667157955e-05 5.865660509574809e-05 -4.6026547951443764e-05;5.865660509574814e-05 0.0005644510772981891 -0.0005648017670874736;-4.602654795144384e-05 -0.0005648017670874737 0.0004531342811589044");
        spin = 1;
    }

    Spin WCHE1
    {
        tensor = matrix("-0.0010009305513396973 0.00020616805963305103 0.00019326373723915018;0.00020616805963305095 -0.00044157311600264787 0.00030726527400815013;0.00019326373723915015 0.0003072652740081502 -0.0003524809407426454");
        spin = 1/2;
    }

    Spin WCHB1
    {
        tensor = matrix("0.0015724517438449976 1.566746046387174e-05 4.7474257073464715e-05;1.566746046387191e-05 0.001515602960922907 6.304611604933337e-05;4.747425707346497e-05 6.304611604933375e-05 0.0017255893265543747");
        spin = 1/2;
    }

	// -------------------------
	// Zeeman interaction
	// -------------------------
	
	Interaction zeeman1
	{
		type = zeeman;
		field = "0 0 0e-5";
		group1 = RPElectron1;
	}	

        Interaction zeeman2
        {
                type = zeeman;
                field = "0.0 0.0 0e-5";
                spins = RPElectron2;
        }

        Interaction zeeman3
        {
                type = zeeman;
                field = "0.0 0.0 0e-5";
                spins = FN5,WCNE1;
		CommonPrefactor = true;
		Prefactor = -0.000109782;
        }

        Interaction zeeman4
        {
                type = zeeman;
                field = "0.0 0.0 0e-5";
                spins = WCHE1,WCHB1;
                CommonPrefactor = true;
                Prefactor = -00151927;
        }

	/ -------------------------
	// Hyperfine interactions
	// -------------------------

    Interaction radical1hyperfine
    {
        type = hyperfine;
        group1 = RPElectron1;
        group2 = FN5;
    }
     Interaction radical2hyperfine
    {
        type = hyperfine;
        group1 = RPElectron2;
        group2 = WCNE1,WCHE1,WCHB1;
    }

	// -------------------------
        // Dipolar interactions
        // -------------------------

        Interaction dipolar
        {
                type = doublespin;
                group1 = RPElectron1;
                group2 = RPElectron2;
		tensor=matrix("3.89e-05 -0.000433 0.000175;-0.000433 -0.000278 0.000251;0.000175 0.000251 0.00024");
		commonprefactor = true;	
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
	Settings general
	{
		steps = 1000;
		notifications = details;
	}

        Action field_strength1
        {
	        type = AddVector;
	        vector = system1.zeeman1.field;
	        direction = "0 0 1";
        	value = 0.01;
        }

        Action field_strength2
        {
	        type = AddVector;
	        vector = system1.zeeman2.field;
	        direction = "0 0 1";
        	value = 0.01;
        }

        Action field_strength3
        {
	        type = AddVector;
	        vector = system1.zeeman3.field;
	        direction = "0 0 1";
        	value = 0.01;
        }

       Action field_strength4
        {
                type = AddVector;
                vector = system1.zeeman4.field;
                direction = "0 0 1";
                value = 0.01;
        }

        // ---------------------------------------------------------
        // Outputs objects
        // ---------------------------------------------------------
        Output orientation
        {
		type = xyz;
                vector = system1.zeeman1.field;
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
		nuclei_list=WCHE1,WCHB1;
	}	
}

