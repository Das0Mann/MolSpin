	    	Interaction Dipolar
	{
	   
		type = DoubleSpin;
        prefactor=0.0020023;

		tensor = matrix("0.039 -0.432 0.175; -0.432 -0.279 0.251; 0.175 0.251 0.239");
		tensortype = ougeneral;
		correlationtime = 100;
		stdev = 0.01;
		modulatedistance=true;
		timestep=1;

        printtensor=true;

		group1 = "electron1";	// Spins in group1 interact with spins in group2
		group2 = "electron2";
	}
    
    Interaction HFI_FN5
	{
	  
		type = DoubleSpin;
        prefactor=0.001;

		tensor = matrix("-0.099 -0.003 0.000; -0.003 -0.087 0.000; 0.000 0.000 1.757");
		tensortype = ougeneral;
		correlationtime = 10;
		modulatedistance=false;

		stdev = 0.01;

        printtensor=false;
	
		timestep=1;
		group1 = "electron1";
		group2 = "FN5";
	}


    Interaction HFI_WNE1
	{
	  
		type = DoubleSpin;
        prefactor=0.001;

		tensor = matrix("-0.053 0.059 -0.046; 0.059  0.564 -0.565; -0.046 -0.565 0.453");
		tensortype = ougeneral;
		correlationtime = 10;
		modulatedistance=false;

		stdev = 0.01;

        printtensor=false;

		timestep=1;
		group1 = "electron1";
		group2 = "WNE1";
	}

