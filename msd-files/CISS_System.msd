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

    //Interaction DipoleDipole
    //{
    //    type = Dipole;
    //    group1 = electron1;
    //    group2 = electron2;
    //    prefactor = 
    //}

    Interaction Zeeman
    {
        type = Zeeman;
        field = "0 0 5e-5"
        spins = electron1,electron2
    }
}