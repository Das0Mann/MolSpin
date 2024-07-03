// -------------------------------------------------------------
// Formic acid formation
// Action vectors and Redfield theory
// -------------------------------------------------------------
SpinSystem system1
{
// ---------------------------------------------------------
// Spins
// ---------------------------------------------------------
Spin RPElectron1
{
        type = electron;
        spin = 1/2;
        tensor = isotropic("2.0007");
}
Spin RPElectron2
{
        type = electron;
        tensor = isotropic("2.0114");
        spin = 1/2;
}
Spin H
{
        spin = 1/2;
        tensor = isotropic("0.05074");
}
// -------------------------
// Zeeman interaction
// -------------------------
Interaction zeeman1
{
        type = zeeman;
        field = "0 0 0.03";                                                                                     // Magnetic field in T
        group1 = RPElectron1;
        tau_c = 0.0089;
	g = 0.0, 0.162057936, 0.162057936, 0.162057936, 0.162057936, 0.162057936;				//Amplitude for spherical tensors (T00, T20)
        //g = 0.0, 0.081000618, 0.081000618, 0.081000618, 0.081000618, 0.081000618;     			// Amplitude for spherical tensors (T00, T20) [For older versions of MolSpin (!), there was an issue with an factor of 2 which made it confusing for the user to get the right amplitude]
        ops = 0;                                                                                                // Spherical tensor basis
        def_g = 1;                                                                                              // Each operator receives seperate amplitude
        def_specdens = 1;                                                                               	// Neglect dynamical shift
        terms = 1;
        coeff = 0;
}

Interaction zeeman2
{
        type = zeeman;
        field = "0.0 0.0 0.03";
        spins = RPElectron2;
}
// -------------------------
// Hyperfine interactions
// -------------------------
Interaction radical1hyperfine
{
        type = hyperfine;
        group1 = RPElectron2;
        group2 = H;
}
// -------------------------
// Spin States
// -------------------------
State Singlet   // |S>
{
        spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;
}
State T0        // |T0>
{
        spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;
}
State Tp        // |T+>
{
        spin(RPElectron2) = |1/2>;
        spin(RPElectron1) = |1/2>;
}
State Tm        // |T->
{
        spin(RPElectron2) = |-1/2>;
        spin(RPElectron1) = |-1/2>;
}
State Identity
{
}
// ---------------------------------------------------------
// Transitions
// ---------------------------------------------------------
Transition Product1
{
        type = sink;
        source = Singlet;
        rate = 0.01;
}
Transition Product3
{
        type = sink;
        source = Singlet;
        rate = 0.00001;
}
Transition Product4
{
        type = sink;
        source = T0;
        rate = 0.00001;
}
Transition Product5
{
        type = sink;
        source = Tp;
        rate = 0.00001;
}
Transition Product6
{
        type = sink;
        source = Tm;
        rate = 0.00001;
}
// -------------------------
// Spin system properties
// -------------------------
Properties properties
{
        initialstate = Identity;
}
}
// -------------------------------------------------------------
Settings
{
Settings general
{
        steps = 88;
        notifications = details;
}
Action field_strength1
{
        type = AddVector;
        vector = system1.zeeman1.field;
        direction = "0 0 1";
        value = 0.01;
}
Action field_strength2
{
        type = AddVector;
        vector = system1.zeeman2.field;
        direction = "0 0 1";
        value = 0.01;
}
Output magnetic_field
{
        type = xyz;
        vector = system1.zeeman1.field;
}

}
// -------------------------------------------------------------
Run
{
Task Method1
{
        type = StaticSS;
        logfile = "example_staticss.log";
        datafile = "example_staticss.dat";
        transitionyields = true;
}
Task Method2
{
        type = Redfield-Relaxation;
        logfile = "redfield-relaxation.log";
        datafile = "redfield-relaxation.dat";
        transitionyields = true;
}
}

