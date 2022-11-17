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
		tensor = matrix("-1.2875202576297838e-05 -2.505399614227668e-06 -5.0497184670114495e-05;-2.5053996142276704e-06 -2.146431872342134e-05 1.1529774968736222e-05;-5.049718467011447e-05 1.1529774968736228e-05 0.0018127464261247248");
		spin = 1;
	}
	Spin FN10
	{
		tensor = matrix("1.4602337753058789e-05 -1.431534329951884e-06 -5.925032405034897e-06;-1.4315343299518834e-06 4.346919744394791e-06 -7.971325515614403e-06;-5.925032405034895e-06 -7.971325515614408e-06 0.00063668607582823");
		spin = 1;
	}
	Spin WHB1
	{
		tensor = matrix("0.00033275634894172845 -2.6301481967838545e-05 3.9318470700047006e-05;-2.6301481967838538e-05 0.00014405219059881547 9.390623041880876e-06;3.9318470700047006e-05 9.390623041880869e-06 0.00015333279574405");
		spin = 1/2;
	}
	Spin WHB2
	{
		tensor = matrix("0.001326552633175887 8.249924579045192e-06 6.38549463292325e-05;8.2499245790452e-06 0.00124974367649421 5.3365442678691767e-05;6.385494632923248e-05 5.336544267869178e-05 0.0013775912985531597");
		spin = 1/2;
	}
//	Spin WNE1
//	{
//		tensor = matrix("7.629450475153465e-05 0.00011391842392384126 -8.881200515586847e-05;0.00011391842392384125 0.0005513061742807311 -0.00044416788274046317;-8.881200515586852e-05 -0.00044416788274046317 0.000445864400672489");
//		spin = 1;
//	}
	Spin WHE1
	{
		tensor = matrix("-0.000634882339148072 0.00024704456358494225 0.00028250282186903486;0.0002470445635849422 -0.00039661775555133226 0.00022053802372068264;0.00028250282186903486 0.00022053802372068264 -0.00018584235691552925");
		spin = 1/2;
	}

	// ---------------------------------------------------------
	// Interactions
	// ---------------------------------------------------------
	Interaction ZeemanElectron
	{
		type = zeeman;
		field = "0 0 0";
		spins = RPElectron1,RPElectron2;
	}

	        Interaction ZeemanNitrogen
        {
                type = zeeman;
                field = "0 0 0";
                spins = FN5,FN10;
        }

                Interaction ZeemanHydrogen
        {
                type = zeeman;
                field = "0 0 0";
                spins = WHB1,WHB2,WHE1;
        }
	

	Interaction radical1hyperfine
	{
		type = hyperfine;
		group1 = RPElectron1;
		group2 = FN5,FN10;
	}
 	Interaction radical2hyperfine
	{
		type = hyperfine;
		group1 = RPElectron2;
		group2 = WHB1,WHB2,WHE1;
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
		source = Singlet;
		rate = 0.001;

	}
	
        Transition Product2
        {
                type = sink;
                source = T0;
                rate = 0.001;

        }

        Transition Product3
        {
                type = sink;
                source = Tp;
                rate = 0.001;

        }

        Transition Product4
        {
                type = sink;
                source = Tm;
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
		steps = 100;
	}
	
	// ---------------------------------------------------------
	// Actions
	// ---------------------------------------------------------

        Action field_strength1
        {
                type = AddVector;
                vector = RPSystem.ZeemanElectron.field;
                direction = "0 0 1";
                value = 0.1;
        }

        Action field_strength2
        {
                type = AddVector;
                vector = RPSystem.ZeemanNitrogen.field;
                direction = "0 0 1";
                value = 0.1;
        }

        Action field_strength3
        {
                type = AddVector;
                vector = RPSystem.ZeemanHydrogen.field;
                direction = "0 0 1";
                value = 0.1;
        }

	
	// ---------------------------------------------------------
	// Outputs objects
	// ---------------------------------------------------------
	Output orientation
	{
		type = vectorxyz;
		vector = RPSystem.ZeemanElectron.field;
	}
}
Run
{
        Task Method1
        {
                type = StaticSS-CIDNP;
                logfile = "example_staticsscidnp.log";
                datafile = "example_staticsscidnp.dat";
                transitionyields = true;
                nuclei_list=N10;
        }
}
