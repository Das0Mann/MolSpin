SpinSystem RPSystem
{
	Spin electron1
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic("2.0023");
	}

	Spin electron2
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic("2.0023");
	}

    Spin N1
	{
		spin = 1/2;///2;
		type = nucleus;
		tensor = isotropic("1.0");
	}
	
	Spin N5
	{
		spin = 1/2;///2;
		type = nucleus;
		tensor = isotropic("1.0");
	}

    Interaction Zeeman
	{
		type = Zeeman;
		field = "5e-02 0 0"; //because timestep is 1e-3
		spins = electron1, electron2;
	}

    Interaction Wc_N5
 	{
		type = Hyperfine;
 		group1 = electron1;
 		group2 = N5;
		prefactor = 6.283185307;//e-03;
		tensor = matrix("-0.36082693 -0.0702137 -1.41518116;
						 -0.0702137  -0.60153649 0.32312139;
      					 -1.41518116  0.32312139 50.80213093");
		CommonPrefactor = false;
		IgnoreTensors = true;
 	}

    Interaction Wc_N1
 	{
 		type = Hyperfine;
 		group1 = electron2;
 		group2 = N1;
		prefactor = 6.283185307;//e-03;
		tensor = matrix(" 2.13814981  3.19255832  -2.48895215;
      					  3.19255832  15.45032887 -12.44778343;
      					 -2.48895215 -12.44778343  12.49532827");
		CommonPrefactor = false;
		IgnoreTensors = true;
 	}

    Interaction DipoleDipole_Wc
	{
		type = Dipole;
		group1 = electron1;
		group2 = electron2;
		tensor = matrix("26.47042689 -55.90357828  50.1679204;
                   		-55.90357828 -20.86385225  76.13493805;
                   		 50.1679204   76.13493805 -5.60657464");
		CommonPrefactor = false;
		IgnoreTensors = true;
	}

    State CISS_Initial
    {
        chi = 0;
        //spins(electron1,electron2) = cos(0.5chi)*(|1/2, -1/2> - |-1/2, 1/2>) + sin(0.5chi)*(|1/2,-1/2> + |-1/2,1/2>);
        spins(electron1,electron2) = 1.414213562|1/2, -1/2>;
    }

    State CISS_Recomb
    {
        chi = 0;
        //spins(electron1,electron2) = cos(0.5chi)*(|1/2, -1/2> - |-1/2, 1/2>) - sin(0.5chi)*(|1/2,-1/2> + |-1/2,1/2>);
        spins(electron1,electron2) = -1.414213562|-1/2, 1/2>;
    }

    State Identity
    {
    }

    Transition spinindependent_decay
    {
        rate = 1;
        source = Identity;
    }

    Transition CISSdecay_Wc
    {
        rate = 1;
        source = CISS_Recomb;
    }

    Properties prop
	{
		initialstate = CISS_Initial;
	}
}

Run
{
    Task CalculateQuantumYeild
    {
        type = StaticSS;
        logfile = "CISS_HoreLauLogTest2.log";
        datafile = "CISS_HoreLauYeildsTest2.dat";
    }
}

Settings
{
    Settings general
    {
        steps = 12500000;
    }

    Action kr
    {
        value = 2;
        type = addscalar;
        scalar = RPSystem.CISSdecay_Wc.rate;
        period = 25000;
    }

    Action kf
    {
        value = 2;
        type = addscalar;
        scalar = RPSystem.spinindependent_decay.rate;
        first = 0;
        last = 25000;
        period = 50;
        loop = true;
    }

    Action scan
    {
        points = 50;
        vector = RPSystem.Zeeman.field;
        type = fibonaccisphere;
        first = 0;
        last = 50;
        loop = true;
    }

    Output field
    {
        type = XYZ;
        vector = RPSystem.Zeeman.field;
    }

    Output krc
    {
        type = scalar;
        scalar = RPSystem.CISSdecay_Wc.rate;
    }

    Output kfc
    {
        type = scalar;
        scalar = RPSystem.spinindependent_decay.rate;
    }
}