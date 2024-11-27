    SpinSystem system1
    {
      // ---------------------------------------------------------
      // Spins
      // ---------------------------------------------------------

        Spin E1
        {
                type = electron;
                spin = 1/2;
                tensor = isotropic(2.0036);                                  // or tensor = isotropic(-2.0023);  g-factor of 2 is default if not specified
        }

        Spin E2
         {
                type = electron;
                spin = 1/2;
                tensor = isotropic(2.0023);
        }

      // -------------------------
      // Zeeman interaction
      // -------------------------
        Interaction zeeman1static
        {
                type = zeeman;
                field = "0.0 0.0 0.337";
                group1 = E1, E2;
        }

        Interaction zeeman2corriolis
        {
                type = zeeman;
                field = "0.0 0.0 -0.3389";
                spins = E1, E2;
        }

        Interaction zeeman3mwstrength
        {
                type = zeeman;
                field = "0.000001 0.0 0.0";
                spins = E1, E2;
        }

         Interaction exchange
         {
                type = doublespin;
                group1 = E1;
                group2 = E2;
                tensor=isotropic("-0.0005");
                ignoretensors=true;
 	        prefactor = 2.0023;
         }

      // -------------------------
      // Spin States
      // -------------------------

        State Singlet	// |S>
        {
                spins(E1,E2) = |1/2,-1/2> - |-1/2,1/2>;
        }

        State T0	// |T0>
        {
                spins(E1,E2) = |1/2,-1/2> + |-1/2,1/2>;
        }
	
	State CISS
	{
      		x = 3.141592654;
		//x = 1.570796327;
		//x = 0.0;
      		spins(E1, E2) = cos(0.5x)(|1/2,-1/2> - |-1/2,1/2>) + sin(0.5x)(|1/2,-1/2> + |-1/2,1/2>);
	}

        State Tp	// |T+>
        {
                spin(E2) = |1/2>;
                spin(E1) = |1/2>;
        }

        State Tm	// |T->
        {
                spin(E2) = |-1/2>;
                spin(E1) = |-1/2>;
        }

        State Identity
        {
        }

        // ---------------------------------------------------------
        // Transitions
        // ---------------------------------------------------------
        
	Transition Product2
        {
                type = sink;
                source = Identity;
                rate = 0.001;
        }

	// Relaxation

	Operator T1RPElectron1 {Type = relaxationt1; rate = 0.032; spins = E1;}
	Operator T2RPElectron1 {Type = relaxationt2; rate = 0.016; spins = E1;}

	Operator T1RPElectron2 {Type = relaxationt1; rate = 0.025; spins = E2;}
	Operator T2RPElectron2 {Type = relaxationt2; rate = 0.012; spins = E2;}


       // -------------------------
       // Spin system properties
       // -------------------------

        Properties properties
        {
                initialstate = CISS;
        }

    }
// -------------------------------------------------------------
    Settings
    {
        Settings general
        {
                steps = 1600;
                notifications = details;
        }

	Action field_strength1 {type = AddVector; vector = system1.zeeman1static.field; direction = "0 0 1"; value = 0.000002;}

        // ---------------------------------------------------------
        // Outputs objects
        // ---------------------------------------------------------

        Output orientation
        {
                type = vectorangle;
                vector = system1.zeeman1static.field;
                reference = "0 0 1";                                    // Get angle relative to z-axis
        }

        Output orientation
        {
                type = xyz;
                vector = system1.zeeman1static.field;
        }

    }

// -------------------------------------------------------------
    Run
    {
        Task Method1
        {
                type = StaticSS-Spectra;
                method = timeinf;                                                        // timeinf or timeevo, timeinf allows to propagate spin density directly to t = inf.
                integration = false;                                         // specify if method = timeevo, if true allows to integrate yields in time from [0 to T] where T lies in [0, totaltime] interval with the step = timestep
                cidsp = false;                                                            // cidsp=true allows to calsulate chemically-induced spin polarization
                spinlist = E1, E2;                                                        // list of the spins to calculate polarization for
                logfile = "test.log";
                datafile = "test.dat";
        }
    }
