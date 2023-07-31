// -------------------------------------------------------------
// FAD/TRP radical pair with motion from simulation of DmCry
// -------------------------------------------------------------
SpinSystem RPSystem
{
// ---------------------------------------------------------
// Spins
// ---------------------------------------------------------
// Based on DmCry (200 ns simulation repeated 50 times)
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

 Spin FADN5
 {
 spin = 1/2;
 }

 // ---------------------------------------------------------
 // Interactions
 // ---------------------------------------------------------


 Interaction radical1hyperfine1
 {
 prefactor = 0.001;
 type = hyperfine;
 group1 = RPElectron1;
 group2 = FADN5;
 tensor = isotropic(5); 
 }

 Interaction zeeman1
 {
 type = zeeman;
 prefactor = 0.001;
 field = "0.0 0.0 0.9" + trajectory("magneticfield.mst");
 spins = RPElectron1, RPElectron2;
 trajectory = "magneticfield.mst";
 }
 
 // ---------------------------------------------------------
 // Spin States
 // ---------------------------------------------------------
 State singlet {spins(RPElectron1,RPElectron2) = |1/2,-1/2> - |-1/2,1/2>;}
 State t0 {spins(RPElectron1,RPElectron2) = |1/2,-1/2> + |-1/2,1/2>;}
 State tp {spin(RPElectron2) = |1/2>;spin(RPElectron1) = |1/2>;}
 State tm {spin(RPElectron2) = |-1/2>;spin(RPElectron1) = |-1/2>;}
 State identity {}
 
 // ---------------------------------------------------------
 // Transitions
 // ---------------------------------------------------------
 Transition Product1
 {
 type = sink;
 source = identity;
 rate = 1e-03;
 }

 // ---------------------------------------------------------
 // SpinSystem Properties
 // ---------------------------------------------------------
 Properties properties
 {
 initialstate = singlet;
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
 time = 0; // Used for the static reference calculation
 }

 // ---------------------------------------------------------
 // Actions
 // ---------------------------------------------------------
 Action scan
 {
 axis = "0 1 0";
 value = 1;
 type = rotatevector;
 vector = RPSystem.zeeman1.field;
 }

 // ---------------------------------------------------------
 // Outputs
 // ---------------------------------------------------------
 // Magnetic field orientation
 Output orientation
 {
 type = vectorangle;
 vector = RPSystem.zeeman1.field;
 reference = "0 0 1";
 }

 // Output objects to check that the trajectory is loaded properly
 Output FADN5_axis1 {type = components;vector = RPSystem.FADN5.Tensor.Axis1;}
 Output FADN5_axis2 {type = components;vector = RPSystem.FADN5.Tensor.Axis2;}
 Output FADN5_axis3 {type = components;vector = RPSystem.FADN5.Tensor.Axis3;}

 Output RPE1_qAxis1 {type = components;vector = RPSystem.RPElectron1.quantizationaxis1;}
 Output RPE1_qAxis2 {type = components;vector = RPSystem.RPElectron1.quantizationaxis2;}
 Output RPE1_qAxis3 {type = components;vector = RPSystem.RPElectron1.quantizationaxis3;}
 }

 Run
 {
 // Calculation using MD trajectory
 Task main
 {
 logfile = "logfile.txt";
 datafile = "results.dat";
 type = DynamicHS-TimeEvolution;
 calculateyields = true;
 timestep = 0.002;
 totaltime = 10; // Trace at total time: 0.67%, i.e. less than 1% left
 }
 }
