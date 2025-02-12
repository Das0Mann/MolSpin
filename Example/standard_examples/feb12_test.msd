```// -------------------------------------------------------------
// MolSpin Input file example
// See the user manual for more information
// Find it at www.molspin.eu
// -------------------------------------------------------------
SpinSystem system1
{
        // Electrons
        Spin E1 {type = electron;spin = 1/2;tensor = isotropic(2.0023);}
        Spin E2 {type = electron;spin = 1/2;tensor = isotropic(2.0023);}

        // Nitrogens (FAD)
        Spin FADN5 {type = nucleus;spin = 1;tensor = isotropic(1.0);}

        // Nitrogens (TRP)
        Spin TRPN1 {type = nucleus;spin = 1;tensor = isotropic(1.0);}

        // -------------------------
        // Zeeman interaction
        // -------------------------

        Interaction Zeeman1 { type = zeeman;field = "0.000 0.000 0.00049";spins = E1;}
        Interaction Zeeman2 { type = zeeman;field = "0.000 0.000 0.00049";spins = E2;}

        // Broadband field
        Interaction zeeman_BB
        {
                type = Zeeman;
                field = "0.0 0.0 0.000000421875";
                spins = E1, E2;

                fieldtype="broadband";
                minfreq=0.00000314;
                maxfreq=0.6283;
                stdev=0.5;
                components=10000;
                randomorientations = false;
        }

        // -------------------------
        // Hyperfine interactions
        // -------------------------

        Interaction FADHYP1 { prefactor = 1.0e-3; type = hyperfine; group1 = E1; group2 = FADN5; tensor = matrix("0.52007444, -0.17396341,  0.77386319; -0.17396341, 0.04412077, -0.24317156; 0.77386319, -0.24317156, 1.14359574");}
        Interaction TRPHYP1 { prefactor = 1.0e-3; type = hyperfine; group1 = E2; group2 = TRPN1; tensor = matrix("0.52785038, 0.35898356, 0.29370566; 0.35898356, 0.29105545, 0.21422574; 0.29370566, 0.21422574, 0.22693535");}

        // -------------------------
        // Dipolar interactions
        // -------------------------
        Interaction Dipolar { prefactor = -2.0023e-3; IgnoreTensors = true; type = doublespin; group1 = E1; group2 = E2; tensor=matrix("0.44628577, 0.38647182, 0.12319905; 0.38647182, -0.13644849, 0.0615624; 0.12319905, 0.0615624, -0.30983728");}

        // ---------------------------------------------------------
        // Spin States
        // ---------------------------------------------------------

        State Singlet {spins(E1,E2) = |1/2,-1/2> - |-1/2,1/2>;}

        State T0 {spins(E1,E2) = |1/2,-1/2> + |-1/2,1/2>;}

        State Tp {spin(E2) = |1/2>; spin(E1) = |1/2>;}

        State Tm {spin(E2) = |-1/2>; spin(E1) = |-1/2>;}

        State Identity {}

        // ---------------------------------------------------------
        // Transitions
        // ---------------------------------------------------------

        Transition Product1 {type = sink; source = Singlet; rate=1.0e-3;}

        Transition Product2 {type = sink; source = T0; rate=1.0e-3;}

        Transition Product3 {type = sink; source = Tp; rate=1.0e-3;}

        Transition Product4 {type = sink; source = Tm; rate=1.0e-3;}

}
// -------------------------------------------------------------
Settings {

        Settings general {steps = 1; notifications = details;}
}
Run
{
        Task Method1
        {
                type = "dynamichs-direct-yields";
                logfile = "smaller_system_BB_static.log";
                datafile = "smaller_system_BB_static.dat";
                initialstate = Singlet;
                seed = 1;
                totaltime = 7000;
                timestep=0.1;
                propagationmethod = "autoexpm";
                precision = "double";
                yieldcorrections = false;
                transitionyields = true;
        }
}