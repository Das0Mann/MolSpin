// -------------------------------------------------------------
// Simple time-evolution calculation on a radical pair.
// Contains the N5 nucleus on FAD. May apply an oscillating magnetic field.
// -------------------------------------------------------------
SpinSystem RPSystem
{

// ---------------------------------------------------------
// Spins
// ---------------------------------------------------------

Spin FADN5
{
	// Hyperfine tensor from:
	// "The quantum needle of the avian magnetic compass."
	// H. G. Hiscock et al. PNAS, 113, 4634-4639 (2016)
tensor = matrix("-0.0989 0.0039 0.0;0.0039 -0.0881 0.0;0.0 0.0 1.7569");
spin = 1;
}

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
tau_c = 1;
g = 1;
ops = 0;
coeff = 1;
}

Interaction zeeman1
{
prefactor = 0.001; // Convert T to mT
type = zeeman;
field = "0.0 0.035355339 0.035355339 ";
spins = RPElectron1,RPElectron2;
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

State Identity // Identity projection
{
}

	Transition spinindependent_decay
	{
		rate = 0;
		source = Identity;
	}

// ---------------------------------------------------------
// SpinSystem Properties
// ---------------------------------------------------------
 Properties properties
 {
 Initialstate = Singlet;
 }
}
// -------------------------------------------------------------
Settings
{
	Settings general
	{
		steps = 1; //88;
		notifications = details;
	}

}
// -------------------------------------------------------------
Run
{
	Task Method1
	{
		type = redfield-relaxation-timeevolution;
                logfile = "redfield-relaxation-timeevolution.log";
                datafile = "redfield-relaxation-timeevolution.dat";
                Timestep = 0.1;
                TotalTime = 1000;
	
	}

        Task Method2
        {
               type = "staticss-timeevolution";
		//type = "StaticHS-Direct-TimeEvo";

                logfile = "example_timeevolution.log";
                datafile = "example_timeevolution.dat";
                Timestep = 0.1;
                TotalTime = 1000;

        }

}

