// -------------------------------------------------------------
SpinSystem RPSystem
{
	// ---------------------------------------------------------
	// Spins
	// ---------------------------------------------------------
	Spin RPElectron1
	{
		type = electron;
		tensor = isotropic(2);
		spin = 1/2;
	}
	
	Spin RPElectron2
	{
		type = electron;
		tensor = isotropic(2);
		spin = 1/2;
	}
	Spin FH1p
	{
		tensor = matrix("0.0002515421957733931 -3.0184136546971e-05 3.053933310435657e-05;-3.0184136546971e-05 0.00013889034867959287 -6.04793366312528e-06;3.0539333104356564e-05 -6.047933663125281e-06 0.0001250599785540463");
		spin = 1/2;
	}
	Spin WHB1
	{
		tensor = matrix("0.00033275634894172845 -2.6301481967838545e-05 3.9318470700047006e-05;-2.6301481967838538e-05 0.00014405219059881547 9.390623041880876e-06;3.9318470700047006e-05 9.390623041880869e-06 0.00015333279574405");
		spin = 1/2;
	}


	// ---------------------------------------------------------
	// Interactions
	// ---------------------------------------------------------
	Interaction zeeman1
	{
		type = zeeman;
		field = "4.357787137e-06 0 4.98097349e-05";
		spins = RPElectron1;
	}

        Interaction zeeman2
        {
                type = zeeman;
                field = "4.357787137e-06 0 4.98097349e-05";
                spins = RPElectron2;
        }


	Interaction radical1hyperfine
	{
		type = hyperfine;
		group1 = RPElectron1;
		group2 = FH1p;
		tau_c = 0;
		g = 0;
	}
 	Interaction radical2hyperfine
	{
		type = hyperfine;
		group1 = RPElectron2;
		group2 = WHB1;
	}
 
	// ---------------------------------------------------------
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
	
	State mp	// |-+>
	{
		spin(RPElectron2) = |1/2>;
		spin(RPElectron1) = |-1/2>;
	}
	
	State pm	// |+->
	{
		spin(RPElectron2) = |-1/2>;
		spin(RPElectron1) = |1/2>;
	}
	
	State Identity	// Identity projection
	{
	}
	
	// ---------------------------------------------------------
	// Transitions
	// ---------------------------------------------------------
	Transition Product1
	{
		type = sink;
		source = Identity;	// spin-independent reaction
		rate = 0.001;

	}
	
	// ---------------------------------------------------------
	// SpinSystem Properties
	// ---------------------------------------------------------
	Properties properties
	{
		Initialstate = Singlet;
	}
}
Settings
{
	// ---------------------------------------------------------
	// General settings
	// ---------------------------------------------------------
	Settings general
	{
		steps = 1;
	}

	// ---------------------------------------------------------
	// Outputs objects
	// ---------------------------------------------------------
	Output orientation
	{
		type = vectorxyz;
		vector = RPSystem.zeeman1.field;
	}
}
Run
{
	Task main
	{
	       type = RP-SymmetricUncoupled-Redfield;	// Specify task class
               logfile = "./5_degrees_logfile.txt";
               datafile = "./5_degrees_result.dat";
	       TimeStep = 1;
	       TotalTime = 16000;
	}
}
