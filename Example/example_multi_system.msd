// -------------------------------------------------------------
// Calculation with two spin systems
// -------------------------------------------------------------
SpinSystem RPSystem1
{
// Spins
Spin FADN5 {tensor = anisotropic("0.087 0.100 1.757");spin = 1;}
Spin TRP {tensor = anisotropic("0 0 1");spin = 1;}
Spin RPElectron1 {type = electron;tensor = isotropic(2);spin = 1/2;}
 Spin RPElectron2 {type = electron;tensor = isotropic(2);spin = 1/2;}

 // Interactions
 Interaction radical1hyperfine {prefactor = 0.001;type = hyperfine;tensor = isotropic(1);group1 = RPElectron1;group2 = FADN5; tau_c = 0.005; g=1; coeff=1;}
 Interaction radical2hyperfine {prefactor = 0.001;type = hyperfine;tensor = isotropic(1);group1 = RPElectron2;group2 = TRP;}
 Interaction zeeman1 {prefactor = 0.001;field = "0.0 0.0 0.05";type = zeeman;spins = RPElectron1,RPElectron2;}

 // Spin States - Radical pair only
 State S {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;}
 State T0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;}
 State Tp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;}
 State Tm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;}
 State Identity {}

 // Spin states - Full specification (RP + FAD + TRP)
 State Spp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State Sp0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State Spm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State S0p {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State S00 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State S0m {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}


 State Smp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State Sm0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State Smm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 State T0pp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State T0p0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State T0pm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State T00p {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State T000 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State T00m {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State T0mp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State T0m0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State T0mm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 State Tppp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State Tpp0 {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State Tppm {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State Tp0p {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State Tp00 {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State Tp0m {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State Tpmp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State Tpm0 {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State Tpmm {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 State Tmpp {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State Tmp0 {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State Tmpm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State Tm0p {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State Tm00 {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State Tm0m {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State Tmmp {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State Tmm0 {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State Tmmm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 // Transitions
 Transition Recombination {type = sink;source = S;rate = 1e-03;}

 #ifdef SINGLE_SYSTEM
 // Products from all sources
 #define SINGLET_PRODUCT_RATE 1e-02
 #define TRIPLET_PRODUCT_RATE 1e-02
 Transition PETSpp {source = Spp;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSp0 {source = Sp0;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSpm {source = Spm;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETS0p {source = S0p;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETS00 {source = S00;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETS0m {source = S0m;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETSmp {source = Smp;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSm0 {source = Sm0;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSmm {source = Smm;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETT0pp {source = T0pp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0p0 {source = T0p0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0pm {source = T0pm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETT00p {source = T00p;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT000 {source = T000;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT00m {source = T00m;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETT0mp {source = T0mp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0m0 {source = T0m0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0mm {source = T0mm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTppp {source = Tppp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTpp0 {source = Tpp0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTppm {source = Tppm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTp0p {source = Tp0p;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTp00 {source = Tp00;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTp0m {source = Tp0m;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTpmp {source = Tpmp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTpm0 {source = Tpm0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTpmm {source = Tpmm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTmpp {source = Tmpp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmp0 {source = Tmp0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmpm {source = Tmpm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTm0p {source = Tm0p;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTm00 {source = Tm00;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTm0m {source = Tm0m;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTmmp {source = Tmmp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmm0 {source = Tmm0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmmm {source = Tmmm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 #else
 #define SINGLET_RATE1 1e-02
 #define TRIPLET_RATE1 1e-02
 // Note: Transitions into mixed ensemble of the TRP nuclei
 Transition ETSpp {source = Spp;targetstate = Spp;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSpp {source = Spp;targetstate = Sp0;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSpp {source = Spp;targetstate = Spm;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETSp0 {source = Sp0;targetstate = Spp;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSp0 {source = Sp0;targetstate = Sp0;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSp0 {source = Sp0;targetstate = Spm;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETSpm {source = Spm;targetstate = Spp;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSpm {source = Spm;targetstate = Sp0;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETSpm {source = Spm;targetstate = Spm;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETS0p {source = S0p;targetstate = S0p;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETS0p {source = S0p;targetstate = S00;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETS0p {source = S0p;targetstate = S0m;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETS00 {source = S00;targetstate = S0p;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETS00 {source = S00;targetstate = S00;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETS00 {source = S00;targetstate = S0m;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETS0m {source = S0m;targetstate = S0p;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETS0m {source = S0m;targetstate = S00;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETS0m {source = S0m;targetstate = S0m;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETSmp {source = Smp;targetstate = Smp;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSmp {source = Smp;targetstate = Sm0;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSmp {source = Smp;targetstate = Smm;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETSm0 {source = Sm0;targetstate = Smp;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSm0 {source = Sm0;targetstate = Sm0;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSm0 {source = Sm0;targetstate = Smm;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETSmm {source = Smm;targetstate = Smp;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSmm {source = Smm;targetstate = Sm0;targetsystem = RPSystem2;rate = SINGLET_RATE1;}
 Transition ETSmm {source = Smm;targetstate = Smm;targetsystem = RPSystem2;rate = SINGLET_RATE1;}

 Transition ETT0pp {source = T0pp;targetstate = T0pp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0pp {source = T0pp;targetstate = T0p0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0pp {source = T0pp;targetstate = T0pm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT0p0 {source = T0p0;targetstate = T0pp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0p0 {source = T0p0;targetstate = T0p0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0p0 {source = T0p0;targetstate = T0pm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT0pm {source = T0pm;targetstate = T0pp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0pm {source = T0pm;targetstate = T0p0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0pm {source = T0pm;targetstate = T0pm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT00p {source = T00p;targetstate = T00p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT00p {source = T00p;targetstate = T000;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT00p {source = T00p;targetstate = T00m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT000 {source = T000;targetstate = T00p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT000 {source = T000;targetstate = T000;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT000 {source = T000;targetstate = T00m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT00m {source = T00m;targetstate = T00p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT00m {source = T00m;targetstate = T000;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT00m {source = T00m;targetstate = T00m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT0mp {source = T0mp;targetstate = T0mp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0mp {source = T0mp;targetstate = T0m0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0mp {source = T0mp;targetstate = T0mm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT0m0 {source = T0m0;targetstate = T0mp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0m0 {source = T0m0;targetstate = T0m0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0m0 {source = T0m0;targetstate = T0mm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETT0mm {source = T0mm;targetstate = T0mp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0mm {source = T0mm;targetstate = T0m0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETT0mm {source = T0mm;targetstate = T0mm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTppp {source = Tppp;targetstate = Tppp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTppp {source = Tppp;targetstate = Tpp0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTppp {source = Tppp;targetstate = Tppm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTpp0 {source = Tpp0;targetstate = Tppp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpp0 {source = Tpp0;targetstate = Tpp0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpp0 {source = Tpp0;targetstate = Tppm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTppm {source = Tppm;targetstate = Tppp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTppm {source = Tppm;targetstate = Tpp0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTppm {source = Tppm;targetstate = Tppm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTp0p {source = Tp0p;targetstate = Tp0p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTp0p {source = Tp0p;targetstate = Tp00;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTp0p {source = Tp0p;targetstate = Tp0m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTp00 {source = Tp00;targetstate = Tp0p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTp00 {source = Tp00;targetstate = Tp00;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTp00 {source = Tp00;targetstate = Tp0m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTp0m {source = Tp0m;targetstate = Tp0p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTp0m {source = Tp0m;targetstate = Tp00;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTp0m {source = Tp0m;targetstate = Tp0m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTpmp {source = Tpmp;targetstate = Tpmp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpmp {source = Tpmp;targetstate = Tpm0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpmp {source = Tpmp;targetstate = Tpmm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTpm0 {source = Tpm0;targetstate = Tpmp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpm0 {source = Tpm0;targetstate = Tpm0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpm0 {source = Tpm0;targetstate = Tpmm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTpmm {source = Tpmm;targetstate = Tpmp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpmm {source = Tpmm;targetstate = Tpm0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTpmm {source = Tpmm;targetstate = Tpmm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTmpp {source = Tmpp;targetstate = Tmpp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmpp {source = Tmpp;targetstate = Tmp0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmpp {source = Tmpp;targetstate = Tmpm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTmp0 {source = Tmp0;targetstate = Tmpp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmp0 {source = Tmp0;targetstate = Tmp0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmp0 {source = Tmp0;targetstate = Tmpm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTmpm {source = Tmpm;targetstate = Tmpp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmpm {source = Tmpm;targetstate = Tmp0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmpm {source = Tmpm;targetstate = Tmpm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTm0p {source = Tm0p;targetstate = Tm0p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTm0p {source = Tm0p;targetstate = Tm00;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTm0p {source = Tm0p;targetstate = Tm0m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTm00 {source = Tm00;targetstate = Tm0p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTm00 {source = Tm00;targetstate = Tm00;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTm00 {source = Tm00;targetstate = Tm0m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTm0m {source = Tm0m;targetstate = Tm0p;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTm0m {source = Tm0m;targetstate = Tm00;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTm0m {source = Tm0m;targetstate = Tm0m;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTmmp {source = Tmmp;targetstate = Tmmp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmmp {source = Tmmp;targetstate = Tmm0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmmp {source = Tmmp;targetstate = Tmmm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTmm0 {source = Tmm0;targetstate = Tmmp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmm0 {source = Tmm0;targetstate = Tmm0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmm0 {source = Tmm0;targetstate = Tmmm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}

 Transition ETTmmm {source = Tmmm;targetstate = Tmmp;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmmm {source = Tmmm;targetstate = Tmm0;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 Transition ETTmmm {source = Tmmm;targetstate = Tmmm;targetsystem = RPSystem2;rate = TRIPLET_RATE1;}
 #endif

 // SpinSystem Properties
 Properties properties {initialstate = S;}
 }
 // -------------------------------------------------------------
 #ifndef SINGLE_SYSTEM
 SpinSystem RPSystem2
 {
 // Spins
 Spin FADN5 {tensor = anisotropic("0.087 0.100 1.757");spin = 1;}
 Spin TRP {tensor = anisotropic("0 0 1");spin = 1;}
 Spin RPElectron1 {type = electron;tensor = isotropic(2);spin = 1/2;}
 Spin RPElectron2 {type = electron;tensor = isotropic(2);spin = 1/2;}

 // Interactions
 Interaction radical1hyperfine {prefactor = 0.001;type = hyperfine;tensor = isotropic(1);group1 = RPElectron1;group2 = FADN5;tau_c = 0.005; g=1; coeff=1;}
 Interaction radical2hyperfine {prefactor = 0.001;type = hyperfine;tensor = isotropic(1);group1 = RPElectron2;group2 = TRP;}
 Interaction zeeman1 {prefactor = 0.001;field = "0.0 0.0 0.05";type = zeeman;spins = RPElectron1,RPElectron2;}

 // Spin States - Radical pair only
 State S {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;}
 State T0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;}
 State Tp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;}
 State Tm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;}
 State Identity {}

 // Spin states - Full specification (RP + FAD + TRP)
 State Spp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State Sp0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State Spm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State S0p {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State S00 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State S0m {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State Smp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State Sm0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State Smm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 State T0pp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State T0p0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State T0pm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State T00p {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State T000 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State T00m {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State T0mp {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State T0m0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State T0mm {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 State Tppp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State Tpp0 {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State Tppm {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State Tp0p {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State Tp00 {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State Tp0m {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State Tpmp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State Tpm0 {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State Tpmm {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 State Tmpp {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |1>;spin(TRP) = |1>;}
 State Tmp0 {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |1>;spin(TRP) = |0>;}
 State Tmpm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |1>;spin(TRP) = |-1>;}

 State Tm0p {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |0>;spin(TRP) = |1>;}
 State Tm00 {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |0>;spin(TRP) = |0>;}
 State Tm0m {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |0>;spin(TRP) = |-1>;}

 State Tmmp {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |-1>;spin(TRP) = |1>;}
 State Tmm0 {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |-1>;spin(TRP) = |0>;}
 State Tmmm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;spin(FADN5) = |-1>;spin(TRP) = |-1>;}

 // Products from all sources
 #define SINGLET_PRODUCT_RATE 1e-02
 #define TRIPLET_PRODUCT_RATE 1e-02
 Transition PETSpp {source = Spp;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSp0 {source = Sp0;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSpm {source = Spm;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETS0p {source = S0p;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETS00 {source = S00;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETS0m {source = S0m;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETSmp {source = Smp;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSm0 {source = Sm0;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}
 Transition PETSmm {source = Smm;targetstate = PS;targetsystem = Products;rate = SINGLET_PRODUCT_RATE;}

 Transition PETT0pp {source = T0pp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0p0 {source = T0p0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0pm {source = T0pm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETT00p {source = T00p;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT000 {source = T000;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT00m {source = T00m;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETT0mp {source = T0mp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0m0 {source = T0m0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETT0mm {source = T0mm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTppp {source = Tppp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTpp0 {source = Tpp0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTppm {source = Tppm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTp0p {source = Tp0p;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTp00 {source = Tp00;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTp0m {source = Tp0m;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTpmp {source = Tpmp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTpm0 {source = Tpm0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTpmm {source = Tpmm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTmpp {source = Tmpp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmp0 {source = Tmp0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmpm {source = Tmpm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTm0p {source = Tm0p;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTm00 {source = Tm00;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTm0m {source = Tm0m;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 Transition PETTmmp {source = Tmmp;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmm0 {source = Tmm0;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}
 Transition PETTmmm {source = Tmmm;targetstate = PT;targetsystem = Products;rate = TRIPLET_PRODUCT_RATE;}

 #define SINGLET_RATE2 1e-02
 #define TRIPLET_RATE2 1e-02
 // Note: Transitions into mixed ensemble of the TRP nuclei
 Transition ETSpp {source = Spp;targetstate = Spp;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSpp {source = Spp;targetstate = Sp0;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSpp {source = Spp;targetstate = Spm;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETSp0 {source = Sp0;targetstate = Spp;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSp0 {source = Sp0;targetstate = Sp0;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSp0 {source = Sp0;targetstate = Spm;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETSpm {source = Spm;targetstate = Spp;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSpm {source = Spm;targetstate = Sp0;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSpm {source = Spm;targetstate = Spm;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETS0p {source = S0p;targetstate = S0p;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETS0p {source = S0p;targetstate = S00;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETS0p {source = S0p;targetstate = S0m;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETS00 {source = S00;targetstate = S0p;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETS00 {source = S00;targetstate = S00;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETS00 {source = S00;targetstate = S0m;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETS0m {source = S0m;targetstate = S0p;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETS0m {source = S0m;targetstate = S00;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETS0m {source = S0m;targetstate = S0m;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETSmp {source = Smp;targetstate = Smp;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSmp {source = Smp;targetstate = Sm0;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSmp {source = Smp;targetstate = Smm;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETSm0 {source = Sm0;targetstate = Smp;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSm0 {source = Sm0;targetstate = Sm0;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSm0 {source = Sm0;targetstate = Smm;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETSmm {source = Smm;targetstate = Smp;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSmm {source = Smm;targetstate = Sm0;targetsystem = RPSystem1;rate = SINGLET_RATE2;}
 Transition ETSmm {source = Smm;targetstate = Smm;targetsystem = RPSystem1;rate = SINGLET_RATE2;}

 Transition ETT0pp {source = T0pp;targetstate = T0pp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0pp {source = T0pp;targetstate = T0p0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0pp {source = T0pp;targetstate = T0pm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT0p0 {source = T0p0;targetstate = T0pp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0p0 {source = T0p0;targetstate = T0p0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0p0 {source = T0p0;targetstate = T0pm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT0pm {source = T0pm;targetstate = T0pp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0pm {source = T0pm;targetstate = T0p0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0pm {source = T0pm;targetstate = T0pm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT00p {source = T00p;targetstate = T00p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT00p {source = T00p;targetstate = T000;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT00p {source = T00p;targetstate = T00m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT000 {source = T000;targetstate = T00p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT000 {source = T000;targetstate = T000;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT000 {source = T000;targetstate = T00m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT00m {source = T00m;targetstate = T00p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT00m {source = T00m;targetstate = T000;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT00m {source = T00m;targetstate = T00m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT0mp {source = T0mp;targetstate = T0mp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0mp {source = T0mp;targetstate = T0m0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0mp {source = T0mp;targetstate = T0mm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT0m0 {source = T0m0;targetstate = T0mp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0m0 {source = T0m0;targetstate = T0m0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT0m0 {source = T0m0;targetstate = T0mm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETT0mm {source = T0mm;targetstate = T0mp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0mm {source = T0mm;targetstate = T0m0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETT0mm {source = T0mm;targetstate = T0mm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTppp {source = Tppp;targetstate = Tppp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTppp {source = Tppp;targetstate = Tpp0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTppp {source = Tppp;targetstate = Tppm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTpp0 {source = Tpp0;targetstate = Tppp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpp0 {source = Tpp0;targetstate = Tpp0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpp0 {source = Tpp0;targetstate = Tppm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTppm {source = Tppm;targetstate = Tppp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTppm {source = Tppm;targetstate = Tpp0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTppm {source = Tppm;targetstate = Tppm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTp0p {source = Tp0p;targetstate = Tp0p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTp0p {source = Tp0p;targetstate = Tp00;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTp0p {source = Tp0p;targetstate = Tp0m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTp00 {source = Tp00;targetstate = Tp0p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTp00 {source = Tp00;targetstate = Tp00;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTp00 {source = Tp00;targetstate = Tp0m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTp0m {source = Tp0m;targetstate = Tp0p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTp0m {source = Tp0m;targetstate = Tp00;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTp0m {source = Tp0m;targetstate = Tp0m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTpmp {source = Tpmp;targetstate = Tpmp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpmp {source = Tpmp;targetstate = Tpm0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpmp {source = Tpmp;targetstate = Tpmm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTpm0 {source = Tpm0;targetstate = Tpmp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpm0 {source = Tpm0;targetstate = Tpm0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpm0 {source = Tpm0;targetstate = Tpmm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTpmm {source = Tpmm;targetstate = Tpmp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpmm {source = Tpmm;targetstate = Tpm0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTpmm {source = Tpmm;targetstate = Tpmm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTmpp {source = Tmpp;targetstate = Tmpp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmpp {source = Tmpp;targetstate = Tmp0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmpp {source = Tmpp;targetstate = Tmpm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTmp0 {source = Tmp0;targetstate = Tmpp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmp0 {source = Tmp0;targetstate = Tmp0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmp0 {source = Tmp0;targetstate = Tmpm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTmpm {source = Tmpm;targetstate = Tmpp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmpm {source = Tmpm;targetstate = Tmp0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmpm {source = Tmpm;targetstate = Tmpm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTm0p {source = Tm0p;targetstate = Tm0p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTm0p {source = Tm0p;targetstate = Tm00;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTm0p {source = Tm0p;targetstate = Tm0m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTm00 {source = Tm00;targetstate = Tm0p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTm00 {source = Tm00;targetstate = Tm00;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTm00 {source = Tm00;targetstate = Tm0m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTm0m {source = Tm0m;targetstate = Tm0p;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTm0m {source = Tm0m;targetstate = Tm00;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTm0m {source = Tm0m;targetstate = Tm0m;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTmmp {source = Tmmp;targetstate = Tmmp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmmp {source = Tmmp;targetstate = Tmm0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmmp {source = Tmmp;targetstate = Tmmm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTmm0 {source = Tmm0;targetstate = Tmmp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmm0 {source = Tmm0;targetstate = Tmm0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmm0 {source = Tmm0;targetstate = Tmmm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 Transition ETTmmm {source = Tmmm;targetstate = Tmmp;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmmm {source = Tmmm;targetstate = Tmm0;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}
 Transition ETTmmm {source = Tmmm;targetstate = Tmmm;targetsystem = RPSystem1;rate = TRIPLET_RATE2;}

 // SpinSystem Properties
 // Start with an empty system: do not specify a spin state
 //Properties properties {initialstate = S;}
 }
 #endif
 // -------------------------------------------------------------
 SpinSystem Products
 {
 Spin s {spin = 1/2;}
 State PS {spin(s) = |1/2>;}
 State PT {spin(s) = |-1/2>;}
 }
 // -------------------------------------------------------------
 // General settings
 Settings
 {
 Settings general {steps = 1;}
 }

 #ifdef SINGLE_SYSTEM
 // Change the name of the output file
 #define results results_singlesystem
 #endif

 // Run section
 Run
 {
 Task main
 {
 timestep = 1e-2;
 totaltime = 100;
 datafile = "results.dat";
 type = multistaticss-redfield-timeevolution;

 }
 }
