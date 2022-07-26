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
	Spin FN5
	{
	        tensor = matrix("-1.0930288604925018e-05 -2.8121800181853453e-06 -6.962551484609009e-05;-2.8121800181853428e-06 -2.222393056626555e-05 1.1500939984266139e-05;-6.962551484609009e-05 1.1500939984266144e-05 0.001792108139667385");
		spin = 1;
	}
	Spin FN10
	{
		tensor = matrix("1.721161329914365e-05 -1.6665051491847646e-06 -2.6062881537238545e-06;-1.6665051491847637e-06 6.5802296548430454e-06 -5.549809506106647e-06;-2.6062881537238596e-06 -5.549809506106646e-06 0.0006395377118263967");
		spin = 1;
	}

	// ---------------------------------------------------------
	// Interactions
	// ---------------------------------------------------------
	Interaction zeeman1
	{
		type = zeeman;
               field = "2.867882181755231e-05 0.0 4.09576022144496e-05";
		spins = RPElectron1,RPElectron2;
	}

	Interaction radical1hyperfine
	{
		type = hyperfine;
		group1 = RPElectron1;
		group2 = FN5;
	}

      Interaction radical3hyperfine
        {
                type = hyperfine;
                group1 = RPElectron1;
                group2 = FN10;
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
	
	State Identity	// Identity projection
	{
	}
	
	// ---------------------------------------------------------
	// Transitions
	// ---------------------------------------------------------
	Transition Product1
	{
		type = sink;
		source = Singlet;	// spin-independent reaction
		rate = 0.001;

	}

	Transition Product2
	{
                type = sink;
                source = T0;      // spin-independent reaction
                rate = 0.001;

	}
	

       Transition Product3
        {
                type = sink;
                source = Tm;      // spin-independent reaction
                rate = 0.001;

        }

       Transition Product4
        {
                type = sink;
                source = Tp;      // spin-independent reaction
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
		steps = 72;
	}
	
	// ---------------------------------------------------------
	// Actions
	// ---------------------------------------------------------
	Action scan
	{
		type = rotatevector;
		vector = RPSystem.zeeman1.field;
               axis = "0 0 1";
		value = 5.0;          
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
	Task static
	{

            type = StaticSS;
               logfile = "logfile_static_35_degrees.txt";
               datafile = "results_static_35_degrees.dat";
            transitionyields = true;

	}

//	Task redfield
//	{
//	        type = redfield-relaxation;
//               logfile = "logfile_redfield_35_degrees.txt";
//               datafile = "results_redfield_35_degrees.dat";
//                transitionyields = true;
//	}
}
