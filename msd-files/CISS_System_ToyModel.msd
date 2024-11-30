SpinSystem RPSystem
{
    Spin electron1
    {
        spin = 1/2;
        type = electron;
        tensor = isotropic(2);
    }

    Spin electron 2
    {
        spin = 1/2;
        type = electron;
        tensor = isotropic(2);
    }

    Spin nucleus1
    {
        spin = 1/2;
        type = nucleus;
        tensor = anisotropic("0.0 0.0 1.5");
    }

    Interaction Hyperfine
    {
        type = Hyperfine;
        group1 = electron1;
        group2 = nucleus1;
        prefactor = 1e-3;
    }

    Interaction DipoleDipole
    {
        type = Dipole;
        group1 = electron1;
        group2 = electron2;
        tensor = isotropic("-0.4e-3")
        ignoretensor = true;
        prefactor = 2.0023;
    }

    Interaction Zeeman
    {
        type = Zeeman;
        field = "0 0 5e-5"
        spins = electron1,electron2
    }

    State Singlet
	{
		spins(electron1, electron2) = |1/2,-1/2> - |-1/2,1/2>;
	}

	State T0
	{
		spins(electron1, electron2) = |1/2,-1/2> + |-1/2,1/2>;
	}

    State TPlus
	{
		spins(electron1) = |1/2>;
		spins(electron2) = |1/2>;
	}

	State TMinus
	{
		spins(electron1) = |-1/2>;
		spins(electron2) = |-1/2>;
	}

    State CISSInitial
    {
        chi = 0;
        spins(electron1,electron2) = cos(0.5chi)*(|1/2, -1/2> - |-1/2, 1/2>) + sin(0.5chi)*(|1/2,-1/2> + |-1/2,1/2>);
    }

    State CISSRecombination
    {
        chi = 0;
        spins(electron1,electron2) = cos(0.5chi)*(|1/2, -1/2> - |-1/2, 1/2>) - sin(0.5chi)*(|1/2,-1/2> + |-1/2,1/2>);
    }

    Transition CISSRecomb
    {
        rate = 1;
        source = CISSRecombination;
    }

    Transition spinindependent_decay
    {
        rate = 1;
        source = Identity;
    }

    Properties prop
    {
        initialstate = CISSInitial;
    }
}

Run
{
    Task CalculateQuantumYeild
    {
        type = StaticSS-TimeEvolution;
        logfile = "CISSToyModelLogFile.log";
        datafile = "CISSToyModelResults.dat";

        timestep = 1;
        totaltime = 12000;
    }
}

Settings
{
    Settings general
	{
		steps = 20;
	}

pi = 3.141592654;
    Action CISS
	{
		value = 0.025pi;
		type = addscalar;
		scalar = RPSystem.CISSInitial.chi;
	}

    Action CISS2
    {
        value = 0.025pi;
        type = addscalar;
        scalar = RPSystem.CISSRecomb.chi
    }

}