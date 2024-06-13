SpinSystem RPSystem
{
//definition of spins
	Spin Nitrogen{ 
		type = nucleus;		
		spin = 1;
		tensor = isotropic(1);}
	
	Spin Prot2{
		type = nucleus;
		spin = 1/2;}
		
	Spin Prot3{
		type = nucleus;
		spin = 1/2;}
		
	Spin RPElectron1{ 
		type = electron; 
		tensor = isotropic(2.0023); 
		spin = 1/2;} // 
		
	Spin RPElectron2{ 
		type = electron; 
		tensor = isotropic(2.0023); 
		spin = 1/2;}


	// Hyperfine
	Interaction radical1hyperfine { 
		#IgnoreTensors = true; 
		#prefactor = 2.0023;
		type = hyperfine;
        	group1 = RPElectron1; 
        	group2 = Nitrogen;
        	//tensor = isotropic(0.000926);
        	tensor = isotropic(0.000926)+anisotropic("-0.0005033 -0.0005033 0.00100667");
        	tau_c = 0.06;
        	g = 0.0, 0.04849, 0.04849, 0.04849, 0.04849, 0.04849;
        	ops = 0; coeff = 0; def_g = 1; terms=1;
        	}
    
        	
        Interaction protonhyperfine {
        	Ignoretensors = true; 
        	Prefactor = 2.0023;
        	type = hyperfine;
        	group1 = RPElectron2; 
        	group2 = Prot2;
        	tensor = isotropic(0.00018);
        	}
        
        Interaction protonhyperfine2 {
        	Ignoretensors = true; 
        	Prefactor=2.0023;
        	type = hyperfine;
        	group1 = RPElectron2; 
        	group2 = Prot3;
        	tensor = isotropic(0.00018);
        	}			
	
	// Zeeman with g anisotropy
	Interaction zeeman1 { 
		type = zeeman; 
		field = "0 0 0.02";
        	spins = RPElectron1;
        	//tau_c = 0.65;
        	//g = 0.0, 0.1467927, 0.1467927, 0.1467927, 0.1467927, 0.1467927;
        	//ops = 0; coeff = 0; def_g = 1; terms=1;
        	}
		
		
	Interaction zeeman2 { 
		type = zeeman;
        	field = "0 0 0.02";
        	spins = RPElectron2;
        	}
	
	
	Interaction exchange { 
		ignoreTensors = true;
        	prefactor = 2.0023;
        	type = doublespin;
        	tensor = isotropic(-0.02);
        	group1 = RPElectron1;
        	group2 = RPElectron2;
        	}

// Spin States for RP only
	State Singlet{spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;}
	State T0 { spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;}
	State Tp { spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;}
	State Tm { spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;}
	State Identity{}

        // ---------------------------------------------------------
        // Transitions
        // ---------------------------------------------------------
        Transition Product1
        {
                type = sink;
                source = Singlet;
                rate = 0;//5e-4 entspricht 5e7 1/s
        }

       Transition Product2
        {
                type = sink;
                source = T0;
                rate = 0; //2e-01 enspricht 1e8 1/s
        }

       Transition Product3
        {
                type = sink;
                source = Tp;
                rate = 0;
        }

       Transition Product4
        {
                type = sink;
                source = Tm;
                rate = 0;
        }


// Additional Relaxation operators
 	//Operator kstd { type = RelaxationDephasing; spins = RPElectron1, RPElectron2; rate = 1e-1;}

// Properties
	Properties properties {Initialstate = Singlet;}

}

Settings{
	Settings general{steps = 1;}
}

// Run section
Run
{

        Task Method1
        {
                type = redfield-relaxation-timeevolution;
                logfile = "redfield-relaxation-timeevolution.log";
                datafile = "output.txt";
                Timestep = 0.2;
                TotalTime = 500;
		transitionyields=true;

        }
}
