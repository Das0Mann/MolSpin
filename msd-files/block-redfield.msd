// -------------------------------------------------------------
 // Simple time-evolution calculation on a radical pair.
 // Contains the N5 nucleus on FAD.
 // May apply an oscillating magnetic field.
 // -------------------------------------------------------------
 SpinSystem RPSystem
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
    Spin FADN5
    {
        // Hyperfine tensor from:
        // "The quantum needle of the avian magnetic compass."
        // H. G. Hiscock et al. PNAS, 113, 4634-4639 (2016)
        tensor = matrix("-0.0989 0.0039 0.0;0.0039 -0.0881 0.0;0.0 0.0 1.7569");
        spin = 1;
    }
    // ---------------------------------------------------------
    // Interactions
    // ---------------------------------------------------------
    Interaction radical1hyperfine
    {
        prefactor = 0.001; // Convert T to mT
        type = hyperfine;
        group1 = RPElectron1;
        group2 = FADN5;
        tensor = isotropic(1);
        tau_c = 0.000; // Correlation time in ns
        g = 1; // Amplitude of correlation function
        coeff = 1; // Use spatial tensors for relaxation process
    }
    Interaction zeeman1
    {
        prefactor = 0.001; // Convert T to mT
        type = zeeman;
        field = "0.0 0.035355339 0.035355339 ";
        spins = RPElectron1,RPElectron2;
        tensor = isotropic(2);
    }
    // ---------------------------------------------------------
    // Spin States
    // ---------------------------------------------------------
    State Singlet // |S>
    {
        spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;
    }
    State T0 // |T0>
    {
        spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;
    }
    State Tp // |T+>
    {
        spin(RPElectron2) = |1/2>;
        spin(RPElectron1) = |1/2>;
    }
    State Tm // |T->
    {
        spin(RPElectron2) = |-1/2>;
        spin(RPElectron1) = |-1/2>;
    }

    Transition SingletDecay
	{
		rate = 0.5e-3; 
		//rate = 0;
		source = Singlet;
	}

    Transition SpinIndependentDecay
    {
        rate = 1e-3;
        source = Identity;
    }

    State Identity // Identity projection
    {
    }
    // ---------------------------------------------------------
    // SpinSystem Properties
    // ---------------------------------------------------------
    Properties properties
    {
        Initialstate = Singlet;
    }

    Subsystem RP1
    {
        spins = RPElectron1, RPElectron2, FADN5;
        Interactions = radical1hyperfine, zeeman1;

        transitions_out = SingletDecay,SpinIndependentDecay; 
    }
 }
 // -------------------------------------------------------------
 Settings
 {
    Settings general
    {
        steps = 5; //88;
        notifications = details;
    }

    Action scan
	{
		points = 5;
		vector = RPSystem.zeeman1.field;
		type = fibonaccisphere;
	}
 }
 // -------------------------------------------------------------
 Run
 {
    Task Method1
    {
        type = StaticSS-TimeEvolution;
        //type = MultiRadicalPairSS-TimeEvolution;
        logfile = "redfield-relaxation-timeevolution5.log";
        datafile = "redfield-timeevolution3.dat";
        Timestep = 0.01; // in ns
        TotalTime = 1000; // in ns
    }

    Task Method2
    {
        //type = StaticSS-TimeEvolution;
        type = MultiRadicalPairSS-TimeEvolution;
        logfile = "redfield-relaxation-timeevolution5.log";
        datafile = "redfield-timeevolution4.dat";
        Timestep = 0.01; // in ns
        TotalTime = 1000; // in ns
    }
 }