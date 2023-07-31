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
		tensor = isotropic(1);
	}
	
	Spin RPElectron2
	{
		type = electron;
		spin = 1/2;
		tensor = isotropic(1);
	}

        Spin H1
        {
                spin = 1/2;
		tensor = isotropic(1);
        }

	/ -------------------------
	// Hyperfine interactions
	// -------------------------

     	Interaction radical1hyperfinehydrogenH1asymmetric
    	{
        	type = hyperfine;
        	group1 = RPElectron1;
        	group2 = H1;
                CommonPrefactor = false;
		ignoretensors = true;
		tensor = matrix("0.0 0.0 0.0;0.0 0.0 0.0;0.00125 0.0  0.005");
    	}

        Interaction radical1hyperfinehydrogenH1symmetric
        {
                type = hyperfine;
                group1 = RPElectron1;
                group2 = H1;
                CommonPrefactor = false;
                ignoretensors = true;
                tensor = matrix("0.0 0.0 0.00125;0.0 0.0 0.0;0.00125 0.0  0.005");
        }

	// -------------------------
        // Dipolar interactions
        // -------------------------

        Interaction dipolar
        {
                type = doublespin;
                group1 = RPElectron1;
                group2 = RPElectron2;
		tensor=matrix("-0.00023193 0.0 0.0; 0.0 -0.00023193 0.0; 0.0 0.0 0.000463867");
		ignoretensors=true;
		CommonPrefactor = false;	
        }

        Interaction exchange
        {
                type = exchange;
                group1 = RPElectron1;
                group2 = RPElectron2;
                tensor=isotropic(0.00054);
                CommonPrefactor = false;
		ignoretensors=true
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
		rate = 5.68E-05;	
                //rate = 0.01;
        }

        Transition Product2
        {
                type = sink;
                source = Identity;
		rate = 5.68E-06;
                //rate = 0.001;
        }


        // -------------------------
	// Spin system properties
	// -------------------------

	Properties properties
	{
		initialstate = T0,Tm,Tp;
	}
}
// -------------------------------------------------------------
Settings
{
	Settings general
	{
		steps = 1;
		notifications = details;
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
                logfile = "example_staticsscidnp_paper_values_31.log";
                datafile = "example_staticsscidnp_paper_values_32.dat";
                transitionyields = true;
		nuclei_list = H1;
	}	
}

