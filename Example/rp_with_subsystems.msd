SpinSystem RPSystem
{
	Spin electron1
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic(2);
	}

	Spin electron2
	{
		spin = 1/2;
		type = electron;
		tensor = isotropic(2);
	}

	Spin nucleus1
	{
		spin = 1/2;
		type = nucleus;
		tensor = isotropic("1.0");
	}
	//
	//Spin nucleus2
	//{
	//	spin = 1/2;
	//	type = nucleus;
	//	tensor = anisotropic("0.5 0.5 0");
	//}
	//
	//Spin nucleus3
	//{
	//	spin = 1/2;
	//	type = nucleus;
	//	tensor = anisotropic("0.0 0.5 0.3");
	//}

	Interaction Hyperfine1
 	{
		type = Hyperfine;
 		group1 = electron1;
 		group2 = nucleus1;
 		prefactor = 1e-3;
		subsystem = RP1;
 	}

 	Interaction Hyperfine2
 	{
 		type = Hyperfine;
 		group1 = electron2;
 		group2 = nucleus2;
 		prefactor = 1e-3;
		subsystem = RP1;
 	}

	Interaction Zeeman
	{
		type = Zeeman;
		field = "0 0 5e-5";
		spins = electron1, electron2;
		subsystem = RP1,RP2;
	}

	Interaction Hyperfine3
	{
		type = Hyperfine;
		group1 = electron1;
		group2 = nucleus1;
		prefactor = 1;
		subsystem = RP2;
	}

	Interaction Hyperfine4
	{
		type = Hyperfine;
		group1 = electron2;
		group2 = nucleus3;
		prefactor = 1;
		subsystem = RP2;
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
	
	State Identity
	{
	}

	Transition SingletDecay
	{
		rate = 1;
		source = Singlet;
		sourcesubsystem = RP1;
	}
	
	Transition spinindependent_decay
	{
		rate = 0.1;
		source = Identity;
		sourcesubsystem = RP1;
	}

	Transition spinindependent_decay2
	{
		rate = 1;
		source = Identity;
		sourcesubsystem = RP2;
	}

	Transition RPtransition
	{
		rate = 1e4;
		source = Identity;
		sourcesubsystem = RP1;
		targetsubsystem = RP2;
		//targetstate = Identity;
	}

	Transition RPtransition2
	{
		rate = 1e4;
		source = Identity;
		sourcesubsystem = RP2;
		targetsubsystem = RP1;
	}

	SubSystem RP1 
	{
		spins = electron1, electron2, nucleus1;
		interactions = Hyperfine1, Hyperfine2;
			
		transitions_out = SingletDecay, spinindependent_decay, RPtransition;
		transitions_in = RPtransition2;
	}

	SubSystem RP2
	{
		spins = electron1, electron2, nucleus1;
		interactions = Hyperfine3, Hyperfine4;

		transitions_out = spinindependent_decay2, RPtransition2;
		transitions_in = RPtransition;
	}
	
	Properties prop
	{
		initialstate = Singlet,Zero;
	}	
		
}


Run
{
	Task CalculateQuantumYeild
	{
		type = MultiRadicalPairSS-TimeEvolution;
		logfile = "logfile.log";
		datafile = "result.dat";
		
		timestep = 1e-4;
		totaltime = 20;

	}
}

Settings
{
	Settings general
	{
		steps = 500;
		timestep = 1e-4;
	}

	Action scan
	{
		value = 0.001;
		direction = "0 0 1";
		type = addvector;
		vector = RPSystem.Zeeman.field;
	}

	Output fieldstrength
	{
		type = length;
		vector = RPSystem.Zeeman.field;
	}
}
