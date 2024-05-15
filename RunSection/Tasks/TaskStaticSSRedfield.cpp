//////////////////////////////////////////////////////////////////////////////
// MolSpin - Redfield Theory Task - developed by Luca Gerhards
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <omp.h>
#include <memory>
#include "TaskStaticSSRedfield.h"
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "Spin.h"
#include "Interaction.h"
#include "ObjectParser.h"
#include "Operator.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskStaticSSRedfield Constructors and Destructor
	// -----------------------------------------------------
	TaskStaticSSRedfield::TaskStaticSSRedfield(const MSDParser::ObjectParser &_parser, const RunSection &_runsection) : BasicTask(_parser, _runsection), reactionOperators(SpinAPI::ReactionOperatorType::Haberkorn), productYieldsOnly(false)
	{
	}

	TaskStaticSSRedfield::~TaskStaticSSRedfield()
	{
	}
	// -----------------------------------------------------
	// TaskStaticSSRedfield protected methods
	// -----------------------------------------------------
	bool TaskStaticSSRedfield::RunLocal()
	{
#pragma omp single
		{
			this->Log() << "Running method StaticSSRedfield." << std::endl;

			// If this is the first step, write first part of header to the data file
			if (this->RunSettings()->CurrentStep() == 1)
			{
				this->WriteHeader(this->Data());
			}
		}
		// ----------------------------------------------------------------
		// SETTING UP SPIN SYSTEM AND GETTING DENSITY OPERATOR
		// ----------------------------------------------------------------
		// Loop through all SpinSystems
		auto systems = this->SpinSystems();
		for (auto i = systems.cbegin(); i < systems.cend(); i++)
		{
			// Make sure we have an initial state
			auto initial_states = (*i)->InitialState();
			if (initial_states.size() < 1)
			{
				this->Log() << "Skipping SpinSystem \"" << (*i)->Name() << "\" as no initial state was specified." << std::endl;
				continue;
			}

			this->Log() << "\nStarting with SpinSystem \"" << (*i)->Name() << "\"." << std::endl;

			// Obtain a SpinSpace to describe the system
			SpinAPI::SpinSpace space(*(*i));
			arma::cx_mat rho0;
			space.UseSuperoperatorSpace(false);
			space.SetReactionOperatorType(this->reactionOperators);

			// Get the initial state
			for (auto j = initial_states.cbegin(); j < initial_states.cend(); j++)
			{
				arma::cx_mat tmp_rho0;
				if (!space.GetState(*j, tmp_rho0))
				{
					this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
					continue;
				}
				if (j == initial_states.cbegin())
					rho0 = tmp_rho0;
				else
					rho0 += tmp_rho0;
			}

			rho0 /= arma::trace(rho0); // The density operator should have a trace of 1

			// ----------------------------------------------------------------
			// GET THE HAMILTONIAN
			// ----------------------------------------------------------------
			arma::cx_mat H;
			if (!space.Hamiltonian(H))
			{
				this->Log() << "Failed to obtain Hamiltonian in superspace." << std::endl;
				continue;
			}

			// ----------------------------------------------------------------
			// DIAGONALIZATION OF H0// We need all of these operators
			// ----------------------------------------------------------------
			arma::cx_mat eigen_vec; // To hold eigenvectors
			arma::vec eigen_val;	// To hold eigenvalues
			arma::cx_mat eig_val_mat;

			this->Log() << "Starting diagonalization..." << std::endl;
			arma::eig_sym(eigen_val, eigen_vec, H);
			this->Log() << "Diagonalization done! Eigenvalues: " << eigen_val.n_elem << ", eigenvectors: " << eigen_vec.n_cols << std::endl;

			// ----------------------------------------------------------------
			// CONSTRUCTING TRANSITION MATRIX "domega" OUT OF EIGENVALUES OF H0
			// ----------------------------------------------------------------
			arma::cx_mat domega(size(H));

			for (int k = 0; k < (int)(eigen_val.size()); k++)
			{
				for (int s = 0; s < (int)(eigen_val.size()); s++)
				{
					domega(k, s) = eigen_val(k) - eigen_val(s);
				}
			}

			// Constructing diagonal matrix with eigenvalues of H0
			eig_val_mat = diagmat(arma::conv_to<arma::cx_mat>::from(eigen_val));

			// Rotate density operator in eigenbasis of H0
			rho0 = (eigen_vec.t() * rho0 * eigen_vec);
			rho0 = rho0 / trace(rho0);

			// ---------------------------------------------------------------
			// SETUP RELAXATION OPERATOR
			// ---------------------------------------------------------------

			// Redfield tensor
			arma::cx_mat R;
			R.set_size(size(kron(H, H)));
			R.zeros();

			// Temporary Redfield tensor
			arma::cx_mat tmp_R;
			tmp_R.set_size(size(kron(H, H)));
			tmp_R.zeros();

			// R Tensor array pointer for parallelization
			arma::cx_mat **ptr_R = NULL;
			int threads = 1;

// Get number of threads for tensor pointer arrays
#pragma omp parallel
			{
				threads = omp_get_num_threads();
			}

			this->Log() << "Setting up threads. " << threads << " CPU are used." << std::endl;

			ptr_R = new arma::cx_mat *[threads];

			for (int l = 0; l < threads; l++)
			{
				ptr_R[l] = new arma::cx_mat(tmp_R);
				*ptr_R[l] *= 0.0;
			}

			// Spectral density matrix
			arma::cx_mat SpecDens;
			SpecDens.set_size(size(domega));
			SpecDens.zeros();

			// Unit matrix
			arma::cx_mat one;
			one.set_size(size(domega));
			one.eye();

			// Correlation function set ups
			std::vector<double> tau_c_list;
			std::vector<double> ampl_list;
			arma::cx_double ampl_combined = 0.0;

			// Defining variables for parameter of user
			int terms = 0;
			int def_g = 0;
			int def_specdens = 0;
			int def_multexpo = 0;
			int ops = 0;
			int coeff = 0;
			int slip = 0;

			this->Log() << "Contructing all required pre-matrices." << std::endl;

			// Defining varibales for loops and storage of operators
			int num_op = 0;
			int num_element = 0;
			arma::cx_mat **ptr_Tensors = NULL;

			// Spin-Operators
			arma::cx_mat *Sz1 = new arma::cx_mat;
			arma::cx_mat *Sx1 = new arma::cx_mat;
			arma::cx_mat *Sy1 = new arma::cx_mat;
			arma::cx_mat *Sz2 = new arma::cx_mat;
			arma::cx_mat *Sx2 = new arma::cx_mat;
			arma::cx_mat *Sy2 = new arma::cx_mat;

			// Double-Spin operators
			arma::cx_mat *Sx1Sx2 = new arma::cx_mat;
			arma::cx_mat *Sx1Sy2 = new arma::cx_mat;
			arma::cx_mat *Sx1Sz2 = new arma::cx_mat;
			arma::cx_mat *Sy1Sx2 = new arma::cx_mat;
			arma::cx_mat *Sy1Sy2 = new arma::cx_mat;
			arma::cx_mat *Sy1Sz2 = new arma::cx_mat;
			arma::cx_mat *Sz1Sx2 = new arma::cx_mat;
			arma::cx_mat *Sz1Sy2 = new arma::cx_mat;
			arma::cx_mat *Sz1Sz2 = new arma::cx_mat;

			// T0 for rank 0 & 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - space.Rk1SphericalTensorXXX)
			arma::cx_mat *T0_rank_0 = new arma::cx_mat;
			arma::cx_mat *T0_rank_2 = new arma::cx_mat;

			// Tp1 & T1m for rank 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - space.Rk1SphericalTensorXXX)
			arma::cx_mat *Tp1 = new arma::cx_mat;
			arma::cx_mat *Tm1 = new arma::cx_mat;

			// Tp2 & Tp2 for rank 2
			arma::cx_mat *Tp2 = new arma::cx_mat;
			arma::cx_mat *Tm2 = new arma::cx_mat;

			// arma::cx_mat ** Tensors_rotating = NULL;
			int k;
			int s;
			int l;

			// ------------------------------------------------------------------
			// STARTING WITH RELAXATION MATRIX CONSTRUCTION - GOOD LUCK
			// ------------------------------------------------------------------

			this->Log() << "Starting with construction of relaxation matrix." << std::endl;
			for (auto interaction = (*i)->interactions_cbegin(); interaction < (*i)->interactions_cend(); interaction++)
			{
				// Chosen parameter in input file
				this->Log() << "------------------------------------------" << std::endl;
				this->Log() << "Chosen input parameter for interaction:  " << (*interaction)->Name() << std::endl;
				this->Log() << "------------------------------------------" << std::endl;
				(*interaction)->Properties()->Get("terms", terms);
				this->Log() << "terms == " << terms << std::endl;
				(*interaction)->Properties()->Get("def_g", def_g);
				this->Log() << "def_g == " << def_g << std::endl;
				(*interaction)->Properties()->Get("def_multexpo", def_multexpo);
				this->Log() << "def_multexpo == " << def_multexpo << std::endl;
				(*interaction)->Properties()->Get("ops", ops);
				this->Log() << "ops == " << ops << std::endl;
				(*interaction)->Properties()->Get("coeff", coeff);
				this->Log() << "coeff == " << coeff << std::endl;
				(*interaction)->Properties()->Get("slip", slip);
				this->Log() << "slip == " << slip << std::endl;
				(*interaction)->Properties()->Get("def_specdens", def_specdens);
				this->Log() << "def_specdens == " << def_specdens << std::endl;
				this->Log() << "------------------------------------------" << std::endl;

				// Check if def_multexpo keyword is used
				if ((*interaction)->Properties()->Get("def_multexpo", def_multexpo) && def_multexpo == 1)
				{
					this->Log() << "Multiple exponential fits found for different operators" << std::endl;
					arma::cx_mat **ptr_SpecDens = NULL;
					arma::mat tau_c_mat;
					arma::mat ampl_mat;

					// Check if tau_c and g lists are available and get the lists in matrix form (tau_c_mat(m,n) - m -> different correlaton function, n -> element in correlation function)
					if ((*interaction)->Properties()->GetMatrix("tau_c", tau_c_mat))
					{
						if ((*interaction)->Properties()->GetMatrix("g", ampl_mat))
						{
							// Check if SingleSpin or DoubleSpin interaction
							if ((*interaction)->Type() == SpinAPI::InteractionType::SingleSpin)
							{
								// ------------------------------------------------------------------
								// Relaxation Matrix Construction for Interactions of Type::SingleSpin
								// ------------------------------------------------------------------

								// Get groups with respect to interaction
								auto group1 = (*interaction)->Group1();
								// Loop through groups to get all interaction type
								for (auto s1 = group1.cbegin(); s1 < group1.cend(); s1++)
								{
									// Which operator basis was chosen:
									if ((*interaction)->Properties()->Get("ops", ops) && ops == 1)
									{
										// --------------------------------------------------------
										// CREATION OF SPIN OPERATORS (Sz, Sx, Sy)
										// --------------------------------------------------------

										this->Log() << "Sx, Sy and Sz operator basis was chosen - ops == 1" << std::endl;

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sz()), *s1, *Sz1))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sx()), *s1, *Sx1))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sy()), *s1, *Sy1))
										{
											return false;
										}

										// Build double-spin operators
										*Sx1Sx2 = (*Sx1); // * Bx
										*Sx1Sy2 = (*Sx1); // * By
										*Sx1Sz2 = (*Sx1); //  *Bz
										*Sy1Sx2 = (*Sy1); // * Bx
										*Sy1Sy2 = (*Sy1); // * By
										*Sy1Sz2 = (*Sy1); // * Bz
										*Sz1Sx2 = (*Sz1); // * Bx
										*Sz1Sy2 = (*Sz1); // * By
										*Sz1Sz2 = (*Sz1); // * Bz

										// Put all tensors on pointer array
										num_op = 9;
										delete[] ptr_Tensors;
										ptr_Tensors = new arma::cx_mat *[num_op];

										ptr_Tensors[0] = Sx1Sx2;
										ptr_Tensors[1] = Sx1Sy2;
										ptr_Tensors[2] = Sx1Sz2;
										ptr_Tensors[3] = Sy1Sx2;
										ptr_Tensors[4] = Sy1Sy2;
										ptr_Tensors[5] = Sy1Sz2;
										ptr_Tensors[6] = Sz1Sx2;
										ptr_Tensors[7] = Sz1Sy2;
										ptr_Tensors[8] = Sz1Sz2;

										// Bx,Bx (Sx1); Bx,By (Sx1); Bx,Bz (Sx1);By, Bz (Sx1);...; Bx,Bx (Sy1)
									}
									else
									{
										// ----------------------------------------------------------------
										// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
										// ----------------------------------------------------------------

										this->Log() << "Irreducible spherical tensor operator basis (Rank 0 & Rank 2) was chosen - ops == 0" << std::endl;

										this->Log() << "Using magentic field to construct SingleSpin irreducible tensors." << std::endl;

										arma::vec static_field = (*interaction)->Field();
										arma::cx_vec complex_field;

										// Make field complex
										complex_field = arma::conv_to<arma::cx_vec>::from(static_field);

										// Rank 0 tensor with m=0
										if (!space.LRk0TensorT0((*s1), complex_field, *T0_rank_0))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and Field"
														<< "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=0
										if (!space.LRk2SphericalTensorT0((*s1), complex_field, *T0_rank_2))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and Field"
														<< "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=1
										if (!space.LRk2SphericalTensorTp1((*s1), complex_field, *Tp1))
										{
											this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name() << " and Field"
														<< "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=-1
										if (!space.LRk2SphericalTensorTm1((*s1), complex_field, *Tm1))
										{
											this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name() << " and Field"
														<< "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=2
										if (!space.LRk2SphericalTensorTp2((*s1), complex_field, *Tp2))
										{
											this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and Field"
														<< "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=-2
										if (!space.LRk2SphericalTensorTm2((*s1), complex_field, *Tm2))
										{
											this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and Field"
														<< "! Skipping." << std::endl;
											continue;
										}

										// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
										num_op = 6;
										delete[] ptr_Tensors;
										ptr_Tensors = new arma::cx_mat *[num_op];
										ptr_Tensors[0] = T0_rank_0;
										ptr_Tensors[1] = T0_rank_2;
										ptr_Tensors[2] = Tm1;
										ptr_Tensors[3] = Tp1;
										ptr_Tensors[4] = Tm2;
										ptr_Tensors[5] = Tp2;

										// Construct spatial spherical Tensors
										SpinAPI::Tensor inTensor(0);
										if ((*interaction)->Properties()->Get("tensor", inTensor) && (*interaction)->Properties()->Get("coeff", coeff) && coeff == 1)
										{
											this->Log() << "Producing spatial shperical tensors with coupling tensor..." << std::endl;
											auto ATensor = (*interaction)->CouplingTensor();
											arma::cx_mat A(3, 3);
											arma::cx_vec Am(9);

											A.zeros();
											A = arma::conv_to<arma::cx_mat>::from((ATensor->LabFrame()));

											Am(0) = (1.0 / sqrt(3.0)) * (A(0, 0) + A(1, 1) + A(2, 2));
											Am(1) = (1.0 / sqrt(6.0)) * (3.0 * A(2, 2) - (A(0, 0) + A(1, 1) + A(2, 2)));
											Am(2) = 0.5 * (A(0, 2) + A(2, 0) + ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
											Am(3) = -0.5 * (A(0, 2) + A(2, 0) - ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
											Am(4) = 0.5 * (A(0, 0) - A(1, 1) - ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));
											Am(5) = 0.5 * (A(0, 0) - A(1, 1) + ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));

											// Rank 0
											*ptr_Tensors[0] = Am(0) * (*ptr_Tensors[0]);
											// Rank 2
											*ptr_Tensors[1] = Am(1) * (*ptr_Tensors[1]);
											*ptr_Tensors[2] = -Am(3) * (*ptr_Tensors[2]);
											*ptr_Tensors[3] = -Am(2) * (*ptr_Tensors[3]);
											*ptr_Tensors[4] = Am(4) * (*ptr_Tensors[4]);
											*ptr_Tensors[5] = Am(5) * (*ptr_Tensors[5]);
										}
										else
										{
											// Norm
											*ptr_Tensors[0] = (*ptr_Tensors[0]);
											*ptr_Tensors[1] = (*ptr_Tensors[1]);
											*ptr_Tensors[2] = -(*ptr_Tensors[2]);
											*ptr_Tensors[3] = -(*ptr_Tensors[3]);
											*ptr_Tensors[4] = (*ptr_Tensors[4]);
											*ptr_Tensors[5] = (*ptr_Tensors[5]);
										}
									}

// Rotate tensors in eigenbasis of H0
#pragma omp parallel for firstprivate(ptr_Tensors) shared(eigen_vec) num_threads(threads)
									for (k = 0; k < num_op; k++)
									{
										*ptr_Tensors[k] = (eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec);
									}

									// Number of elments for SpecDens. Important for delete statemant later
									if ((*interaction)->Properties()->Get("terms", terms) && terms == 0)
									{
										num_element = num_op * num_op;
									}
									else if ((*interaction)->Properties()->Get("terms", terms) && terms == 1)
									{
										num_element = num_op;
									}

									// Contructing dimension of ptr_SpecDens
									delete[] ptr_SpecDens;
									ptr_SpecDens = new arma::cx_mat *[num_element];

									for (int l = 0; l < num_element; l++)
									{
										ptr_SpecDens[l] = new arma::cx_mat(domega);
										*ptr_SpecDens[l] *= 0.0;
									}

									// TODO: Slippage Operator for multiexponential loop

									this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << std::endl;

									if ((*interaction)->Properties()->Get("terms", terms) && terms == 0 && ops == 1)
									{
										this->Log() << "Setting up Redfield tensor for each component (Axxxx, Axxxy, Axxxz, Axxyx, ...)" << std::endl;

										// Counter for inner loop
										int m = 0;
										// Counter to check if some rows are 0 and must not be computed
										double max_row_value;
										// Do the autocorrelation terms for each of the 9 terms in the tensor
										for (k = 0; k < num_op; k++)
										{
											for (s = 0; s < num_op; s++)
											{
												max_row_value = max(ampl_mat.row(m));

												// Check if the maximum value is zero
												if (max_row_value == 0.0)
												{
													this->Log() << "The maximum value of row (correlation function) " << m << " is zero." << std::endl;
												}
												else
												{
													for (int n = 0; n < (int)ampl_mat.n_cols; n++)
													{
														// ----------------------------------------------------------------
														// CONSTRUCTING SPECTRAL DENSITY MATRIX
														// ----------------------------------------------------------------

														if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
														{
															SpecDens *= 0.0;

															if (!ConstructSpecDensSpecific(1, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
															{
																this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
																continue;
															}

															SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
														}
														else
														{
															SpecDens *= 0.0;

															if (!ConstructSpecDensSpecific(0, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
															{
																this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
																continue;
															}

															SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
														}
														*ptr_SpecDens[m] += SpecDens;
													}

													// -----------------------------------------------------------------
													// CONSTRUCTING R MATRIX
													// -----------------------------------------------------------------

													tmp_R *= 0.0;

													if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[s]), *ptr_SpecDens[m], tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													R += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;

														if (!Redfieldtensor((*ptr_Tensors[s]), (*ptr_Tensors[k]), *ptr_SpecDens[m], tmp_R))
														{
															this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
															continue;
														}

														R += tmp_R;
													}
												}

												m = m + 1;
											}
										}
									}
									else if ((*interaction)->Properties()->Get("terms", terms) && terms == 1 && ops == 1)
									{
										this->Log() << "Setting up Redfield tensor just for each operator seperartely. No cross-correlation." << std::endl;

										// Counter for inner loop
										int m = 0;
										// Counter to check if some rows are 0 and must not be computed
										double max_row_value;

										// Do the autocorrelation terms for each of the 3 terms
										for (k = 0; k < num_op; k++)
										{
											max_row_value = max(ampl_mat.row(m));

											// Check if the maximum value is zero
											if (max_row_value == 0.0)
											{
												this->Log() << "The maximum value of row (correlation function) " << m << " is zero." << std::endl;
											}
											else
											{
												for (int n = 0; n < (int)ampl_mat.n_cols; n++)
												{
													// ----------------------------------------------------------------
													// CONSTRUCTING SPECTRAL DENSITY MATRIX
													// ----------------------------------------------------------------

													if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
													{
														SpecDens *= 0.0;

														if (!ConstructSpecDensSpecific(1, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}

														SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
													}
													else
													{
														SpecDens *= 0.0;

														if (!ConstructSpecDensSpecific(0, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}

														SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
													}

													*ptr_SpecDens[m] += SpecDens;
												}

												// -----------------------------------------------------------------
												// CONSTRUCTING R MATRIX
												// -----------------------------------------------------------------

												tmp_R *= 0.0;

												if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[k]), *ptr_SpecDens[m], tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												R += tmp_R;
											}

											m = m + 1;
										}
									}
									else
									{
										this->Log() << "Your current options are not implemented yet. Please ensure that you use Cartesian operators (ops = 1)." << std::endl;
									}
								}
							}
							else
							{
								// Get groups with respect to interaction
								auto group1 = (*interaction)->Group1();
								auto group2 = (*interaction)->Group2();

								// Loop through groups to get all interaction [group1 = electrons, group2 = particles electrons interact with, commonly nuclei]
								for (auto s1 = group1.cbegin(); s1 < group1.cend(); s1++)
								{
									for (auto s2 = group2.cbegin(); s2 < group2.cend(); s2++)
									{
										if ((*interaction)->Properties()->Get("ops", ops) && ops == 1)
										{

											// --------------------------------------------------------
											// CREATION OF SPIN OPERATORS (Sz, Sx, Sy)
											// --------------------------------------------------------

											this->Log() << "Sx, Sy and Sz operator basis was chosen - ops == 1" << std::endl;

											if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sz()), *s1, *Sz1))
											{
												return false;
											}

											if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sx()), *s1, *Sx1))
											{
												return false;
											}

											if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sy()), *s1, *Sy1))
											{
												return false;
											}

											if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s2)->Sz()), *s2, *Sz2))
											{
												return false;
											}

											if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s2)->Sx()), *s2, *Sx2))
											{
												return false;
											}

											if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s2)->Sy()), *s2, *Sy2))
											{
												return false;
											}

											// Build double-spin operators
											*Sx1Sx2 = (*Sx1) * (*Sx2);
											*Sx1Sy2 = (*Sx1) * (*Sy2);
											*Sx1Sz2 = (*Sx1) * (*Sz2);
											*Sy1Sx2 = (*Sy1) * (*Sx2);
											*Sy1Sy2 = (*Sy1) * (*Sy2);
											*Sy1Sz2 = (*Sy1) * (*Sz2);
											*Sz1Sx2 = (*Sz1) * (*Sx2);
											*Sz1Sy2 = (*Sz1) * (*Sy2);
											*Sz1Sz2 = (*Sz1) * (*Sz2);

											// Put all tensors on pointer array
											num_op = 9;
											delete[] ptr_Tensors;
											ptr_Tensors = new arma::cx_mat *[num_op];

											ptr_Tensors[0] = Sx1Sx2;
											ptr_Tensors[1] = Sx1Sy2;
											ptr_Tensors[2] = Sx1Sz2;
											ptr_Tensors[3] = Sy1Sx2;
											ptr_Tensors[4] = Sy1Sy2;
											ptr_Tensors[5] = Sy1Sz2;
											ptr_Tensors[6] = Sz1Sx2;
											ptr_Tensors[7] = Sz1Sy2;
											ptr_Tensors[8] = Sz1Sz2;
										}
										else
										{
											// -------------------------------------------------------------------
											// Relaxation Matrix Construction for Interactions of Type::DoubleSpin
											// -------------------------------------------------------------------

											this->Log() << "Irreducible spherical tensor operator basis (Rank 0 & Rank 2) was chosen - ops == 0" << std::endl;

											// ----------------------------------------------------------------
											// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
											// ----------------------------------------------------------------

											if (!space.BlRk0TensorT0(*s1, *s2, *T0_rank_0))
											{
												this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Rank 2 tensor with m=0
											if (!space.BlRk2SphericalTensorT0(*s1, *s2, *T0_rank_2))
											{
												this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Rank 2 tensor with m=1
											if (!space.BlRk2SphericalTensorTp1(*s1, *s2, *Tp1))
											{
												this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Rank 2 tensor with m=-1
											if (!space.BlRk2SphericalTensorTm1(*s1, *s2, *Tm1))
											{
												this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Rank 2 tensor with m=2
											if (!space.BlRk2SphericalTensorTp2(*s1, *s2, *Tp2))
											{
												this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Rank 2 tensor with m=-2
											if (!space.BlRk2SphericalTensorTm2(*s1, *s2, *Tm2))
											{
												this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
											num_op = 6;
											delete[] ptr_Tensors;
											ptr_Tensors = new arma::cx_mat *[num_op];
											ptr_Tensors[0] = T0_rank_0;
											ptr_Tensors[1] = T0_rank_2;
											ptr_Tensors[2] = Tm1;
											ptr_Tensors[3] = Tp1;
											ptr_Tensors[4] = Tm2;
											ptr_Tensors[5] = Tp2;

											// Construct spatial spherical Tensors
											SpinAPI::Tensor inTensor(0);
											if ((*interaction)->Properties()->Get("tensor", inTensor) && (*interaction)->Properties()->Get("coeff", coeff) && coeff == 1)
											{
												this->Log() << "Producing spatial shperical tensors with coupling tensor..." << std::endl;
												auto ATensor = (*interaction)->CouplingTensor();
												arma::cx_mat A(3, 3);
												arma::cx_vec Am(6);

												A.zeros();
												A = arma::conv_to<arma::cx_mat>::from((ATensor->LabFrame()));

												Am(0) = (1.0 / sqrt(3.0)) * (A(0, 0) + A(1, 1) + A(2, 2));
												Am(1) = (1.0 / sqrt(6.0)) * (3.0 * A(2, 2) - (A(0, 0) + A(1, 1) + A(2, 2)));
												Am(2) = 0.5 * (A(0, 2) + A(2, 0) - ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
												Am(3) = -0.5 * (A(0, 2) + A(2, 0) + ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
												Am(4) = 0.5 * (A(0, 0) - A(1, 1) - ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));
												Am(5) = 0.5 * (A(0, 0) - A(1, 1) + ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));

												// Rank 0
												*ptr_Tensors[0] = Am(0) * (*ptr_Tensors[0]);
												// Rank 2
												*ptr_Tensors[1] = Am(1) * (*ptr_Tensors[1]);
												*ptr_Tensors[2] = -Am(3) * (*ptr_Tensors[2]);
												*ptr_Tensors[3] = -Am(2) * (*ptr_Tensors[3]);
												*ptr_Tensors[4] = Am(4) * (*ptr_Tensors[4]);
												*ptr_Tensors[5] = Am(5) * (*ptr_Tensors[5]);
											}
											else
											{
												// Norm
												*ptr_Tensors[0] = (*ptr_Tensors[0]);
												*ptr_Tensors[1] = (*ptr_Tensors[1]);
												*ptr_Tensors[2] = -(*ptr_Tensors[2]);
												*ptr_Tensors[3] = -(*ptr_Tensors[3]);
												*ptr_Tensors[4] = (*ptr_Tensors[4]);
												*ptr_Tensors[5] = (*ptr_Tensors[5]);
											}
										}

// Rotate tensors in eigenbasis of H0
#pragma omp parallel for firstprivate(ptr_Tensors) shared(eigen_vec) num_threads(threads)
										for (k = 0; k < num_op; k++)
										{
											*ptr_Tensors[k] = eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec;
										}

										// Number of elments for SpecDens. Important for delete statemant later
										if ((*interaction)->Properties()->Get("terms", terms) && terms == 0)
										{
											num_element = num_op * num_op;
										}
										else if ((*interaction)->Properties()->Get("terms", terms) && terms == 1)
										{
											num_element = num_op;
										}

										// Contructing dimension of ptr_SpecDens
										delete[] ptr_SpecDens;
										ptr_SpecDens = new arma::cx_mat *[num_element];

										for (int l = 0; l < num_element; l++)
										{
											ptr_SpecDens[l] = new arma::cx_mat(domega);
											*ptr_SpecDens[l] *= 0.0;
										}

										// TODO: Slippage Operator for multiexponential loop

										this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << " and spin " << (*s2)->Name() << std::endl;

										if ((*interaction)->Properties()->Get("terms", terms) && terms == 0 && ops == 1)
										{
											this->Log() << "Setting up Redfield tensor for each component (Axxxx, Axxxy, Axxxz, Axxyx, ...)" << std::endl;

											// Counter for inner loop
											int m = 0;
											// Counter to check if some rows are 0 and must not be computed
											double max_row_value;
											// Do the autocorrelation terms for each of the 9 terms in the tensor
											for (k = 0; k < num_op; k++)
											{
												for (s = 0; s < num_op; s++)
												{
													max_row_value = max(ampl_mat.row(m));

													// Check if the maximum value is zero
													if (max_row_value == 0.0)
													{
														this->Log() << "The maximum value of row (correlation function) " << m << " is zero." << std::endl;
													}
													else
													{
														for (int n = 0; n < (int)ampl_mat.n_cols; n++)
														{
															// ----------------------------------------------------------------
															// CONSTRUCTING SPECTRAL DENSITY MATRIX
															// ----------------------------------------------------------------

															if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
															{
																SpecDens *= 0.0;

																if (!ConstructSpecDensSpecific(1, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
																{
																	this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
																	continue;
																}

																SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());

																for (int kok = 0; kok < (int)SpecDens.n_rows; kok++)
																{
																	for (int kek = 0; kek < (int)SpecDens.n_cols; kek++)
																	{
																		SpecDens(kok, kek) = std::real(SpecDens(kok, kek));
																	}
																}
															}
															else
															{
																SpecDens *= 0.0;

																if (!ConstructSpecDensSpecific(0, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
																{
																	this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
																	continue;
																}

																SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
															}
															*ptr_SpecDens[m] += SpecDens;
														}

														std::cout << *ptr_SpecDens[m] << std::endl;

														// -----------------------------------------------------------------
														// CONSTRUCTING R MATRIX
														// -----------------------------------------------------------------

														tmp_R *= 0.0;

														if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[s]), *ptr_SpecDens[m], tmp_R))
														{
															this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
															continue;
														}

														R += tmp_R;

														if (k != s)
														{
															tmp_R *= 0.0;

															if (!Redfieldtensor((*ptr_Tensors[s]), (*ptr_Tensors[k]), *ptr_SpecDens[m], tmp_R))
															{
																this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
																continue;
															}

															R += tmp_R;
														}
													}

													m = m + 1;
												}
											}
										}
										else if ((*interaction)->Properties()->Get("terms", terms) && terms == 1 && ops == 1)
										{

											this->Log() << "Setting up Redfield tensor just for each operator seperartely. No cross-correlation." << std::endl;

											// Counter for inner loop
											int m = 0;
											// Counter to check if some rows are 0 and must not be computed
											double max_row_value;

											// Do the autocorrelation terms for each of the 3 terms
											for (k = 0; k < num_op; k++)
											{
												max_row_value = max(ampl_mat.row(m));

												// Check if the maximum value is zero
												if (max_row_value == 0.0)
												{
													this->Log() << "The maximum value of row (correlation function) " << m << " is zero." << std::endl;
												}
												else
												{
													for (int n = 0; n < (int)ampl_mat.n_cols; n++)
													{
														// ----------------------------------------------------------------
														// CONSTRUCTING SPECTRAL DENSITY MATRIX
														// ----------------------------------------------------------------

														if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
														{
															SpecDens *= 0.0;

															if (!ConstructSpecDensSpecific(1, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
															{
																this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
																continue;
															}

															SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
														}
														else
														{
															SpecDens *= 0.0;

															if (!ConstructSpecDensSpecific(0, static_cast<std::complex<double>>(ampl_mat(m, n)), static_cast<std::complex<double>>(tau_c_mat(m, n)), domega, SpecDens))
															{
																this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
																continue;
															}

															SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
														}

														*ptr_SpecDens[m] += SpecDens;
													}

													// -----------------------------------------------------------------
													// CONSTRUCTING R MATRIX
													// -----------------------------------------------------------------

													tmp_R *= 0.0;

													if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[k]), *ptr_SpecDens[m], tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													R += tmp_R;
												}

												m = m + 1;
											}
										}
										else
										{
											this->Log() << "Your current options are not implemented yet. Please ensure that you use Cartesian operators (ops = 1)." << std::endl;
										}

										this->Log() << "Added relaxation matrix term for interaction " << (*interaction)->Name() << " between spins " << (*s1)->Name() << " and " << (*s2)->Name() << "." << std::endl;
									}
								}
							}
						}
					}

					if (num_element == num_op * num_op)
					{
						for (int i = 0; i < num_op * num_op; i++)
						{
							delete ptr_SpecDens[i]; // Delete individual arma::cx_mat objects
						}
						delete[] ptr_SpecDens; // Delete the array of arma::cx_mat pointers
					}
					else if (num_element == num_op)
					{
						for (int i = 0; i < num_element; i++)
						{
							delete ptr_SpecDens[i]; // Delete individual arma::cx_mat objects
						}
						delete[] ptr_SpecDens; // Delete the array of arma::cx_mat pointers
					}
					else
					{
						// Handle the case when neither l is 9 nor 3, or provide appropriate error handling
						// ...
					}
				}
				// Normal loops without multiple multiexponential lists
				else if ((*interaction)->Properties()->GetList("tau_c", tau_c_list))
				{
					if ((*interaction)->Properties()->GetList("g", ampl_list))
					{
						if ((*interaction)->Type() == SpinAPI::InteractionType::SingleSpin)
						{
							// ------------------------------------------------------------------
							// Relaxation Matrix Construction for Interactions of Type::SingleSpin
							// ------------------------------------------------------------------

							// Get groups with respect to interaction
							auto group1 = (*interaction)->Group1();

							// Loop through groups to get all interaction type
							for (auto s1 = group1.cbegin(); s1 < group1.cend(); s1++)
							{
								// Which operator basis was chosen:
								if ((*interaction)->Properties()->Get("ops", ops) && ops == 1)
								{
									// --------------------------------------------------------
									// CREATION OF SPIN OPERATORS (Sz, Sx, Sy)
									// --------------------------------------------------------

									this->Log() << "Sx, Sy and Sz operator basis was chosen - ops == 1" << std::endl;

									if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sz()), *s1, *Sz1))
									{
										return false;
									}

									if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sx()), *s1, *Sx1))
									{
										return false;
									}

									if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sy()), *s1, *Sy1))
									{
										return false;
									}

									// Put all tensors on pointer array
									num_op = 3;
									delete[] ptr_Tensors;
									ptr_Tensors = new arma::cx_mat *[num_op];
									ptr_Tensors[0] = Sx1;
									ptr_Tensors[1] = Sy1;
									ptr_Tensors[2] = Sz1;
								}
								else
								{
									// ----------------------------------------------------------------
									// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
									// ----------------------------------------------------------------

									this->Log() << "Irreducible spherical tensor operator basis (Rank 0 & Rank 2) was chosen - ops == 0" << std::endl;

									this->Log() << "Using magentic field to construct SingleSpin irreducible tensors." << std::endl;

									arma::vec static_field = (*interaction)->Field();
									arma::cx_vec complex_field;

									// Make field complex
									complex_field = arma::conv_to<arma::cx_vec>::from(static_field);

									// Rank 0 tensor with m=0
									if (!space.LRk0TensorT0((*s1), complex_field, *T0_rank_0))
									{
										this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and Field"
													<< "! Skipping." << std::endl;
										continue;
									}

									// Rank 2 tensor with m=0
									if (!space.LRk2SphericalTensorT0((*s1), complex_field, *T0_rank_2))
									{
										this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and Field"
													<< "! Skipping." << std::endl;
										continue;
									}

									// Rank 2 tensor with m=1
									if (!space.LRk2SphericalTensorTp1((*s1), complex_field, *Tp1))
									{
										this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name() << " and Field"
													<< "! Skipping." << std::endl;
										continue;
									}

									// Rank 2 tensor with m=-1
									if (!space.LRk2SphericalTensorTm1((*s1), complex_field, *Tm1))
									{
										this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name() << " and Field"
													<< "! Skipping." << std::endl;
										continue;
									}

									// Rank 2 tensor with m=2
									if (!space.LRk2SphericalTensorTp2((*s1), complex_field, *Tp2))
									{
										this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and Field"
													<< "! Skipping." << std::endl;
										continue;
									}

									// Rank 2 tensor with m=-2
									if (!space.LRk2SphericalTensorTm2((*s1), complex_field, *Tm2))
									{
										this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and Field"
													<< "! Skipping." << std::endl;
										continue;
									}

									// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
									num_op = 6;
									delete[] ptr_Tensors;
									ptr_Tensors = new arma::cx_mat *[num_op];
									ptr_Tensors[0] = T0_rank_0;
									ptr_Tensors[1] = T0_rank_2;
									ptr_Tensors[2] = Tm1;
									ptr_Tensors[3] = Tp1;
									ptr_Tensors[4] = Tm2;
									ptr_Tensors[5] = Tp2;

									// Construct spatial spherical Tensors
									SpinAPI::Tensor inTensor(0);
									if ((*interaction)->Properties()->Get("tensor", inTensor) && (*interaction)->Properties()->Get("coeff", coeff) && coeff == 1)
									{
										this->Log() << "Producing spatial shperical tensors with coupling tensor..." << std::endl;
										auto ATensor = (*interaction)->CouplingTensor();
										arma::cx_mat A(3, 3);
										arma::cx_vec Am(9);

										A.zeros();
										A = arma::conv_to<arma::cx_mat>::from((ATensor->LabFrame()));

										Am(0) = (1.0 / sqrt(3.0)) * (A(0, 0) + A(1, 1) + A(2, 2));
										Am(1) = (1.0 / sqrt(6.0)) * (3.0 * A(2, 2) - (A(0, 0) + A(1, 1) + A(2, 2)));
										Am(2) = 0.5 * (A(0, 2) + A(2, 0) + ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
										Am(3) = -0.5 * (A(0, 2) + A(2, 0) - ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
										Am(4) = 0.5 * (A(0, 0) - A(1, 1) - ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));
										Am(5) = 0.5 * (A(0, 0) - A(1, 1) + ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));

										// Rank 0
										*ptr_Tensors[0] = Am(0) * (*ptr_Tensors[0]);
										// Rank 2
										*ptr_Tensors[1] = Am(1) * (*ptr_Tensors[1]);
										*ptr_Tensors[2] = -Am(3) * (*ptr_Tensors[2]);
										*ptr_Tensors[3] = -Am(2) * (*ptr_Tensors[3]);
										*ptr_Tensors[4] = Am(4) * (*ptr_Tensors[4]);
										*ptr_Tensors[5] = Am(5) * (*ptr_Tensors[5]);
									}
									else
									{
										// Norm
										*ptr_Tensors[0] = (*ptr_Tensors[0]);
										*ptr_Tensors[1] = (*ptr_Tensors[1]);
										*ptr_Tensors[2] = -(*ptr_Tensors[2]);
										*ptr_Tensors[3] = -(*ptr_Tensors[3]);
										*ptr_Tensors[4] = (*ptr_Tensors[4]);
										*ptr_Tensors[5] = (*ptr_Tensors[5]);
									}
								}

// Rotate tensors in eigenbasis of H0
#pragma omp parallel for firstprivate(ptr_Tensors) shared(eigen_vec) num_threads(threads)
								for (k = 0; k < num_op; k++)
								{
									*ptr_Tensors[k] = (eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec);
								}

								if ((*interaction)->Properties()->Get("slip", slip) && slip == 1)
								{
									this->Log() << "Calculating slippage of inital condition for interaction" << std::endl;
									arma::cx_mat rho0_new;
									rho0_new.set_size(size(rho0));
									rho0_new.zeros();

									if (!Slippage(ptr_Tensors, num_op, eig_val_mat, domega, rho0, static_cast<std::complex<double>>(tau_c_list[0]), rho0_new))
									{
										this->Log() << "There are problems with the construction of the intial condition for slippage - Please check your input." << std::endl;
										continue;
									}

									rho0 = rho0_new;
									rho0_new *= 0.0;
								}

								this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << std::endl;

								if ((*interaction)->Properties()->Get("terms", terms) && terms == 1)
								{
									this->Log() << "No cross-relaxation terms are requested" << std::endl;

									if ((*interaction)->Properties()->Get("def_g", def_g) && def_g == 1)
									{
										this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;

// loop through all operators and construct R tensor in Liouville space
#pragma omp parallel for num_threads(threads)
										for (int l = 0; l < threads; l++)
										{
											*ptr_R[l] *= 0.0;
										}

#pragma omp parallel for firstprivate(tmp_R, SpecDens, ampl_combined) shared(ampl_list, tau_c_list, domega, num_op, ptr_Tensors, ptr_R, interaction, def_specdens) num_threads(threads)
										for (k = 0; k < num_op; k++)
										{
											// ----------------------------------------------------------------
											// CONSTRUCTING SPECTRAL DENSITY MATRIX
											// ----------------------------------------------------------------
											ampl_combined = ampl_list[k] * ampl_list[k];

											if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
											{
												SpecDens *= 0.0;

												if (!ConstructSpecDensSpecific(1, static_cast<std::complex<double>>(ampl_combined), static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
											}
											else
											{
												SpecDens *= 0.0;

												if (!ConstructSpecDensSpecific(0, ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
											}

											// -----------------------------------------------------------------
											// CONSTRUCTING R MATRIX
											// -----------------------------------------------------------------

											tmp_R *= 0.0;

											if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[k]), SpecDens, tmp_R))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											*ptr_R[omp_get_thread_num()] += tmp_R;
										}
									}
									else
									{
										this->Log() << "J is generally constructed for all operators - def_g == 0 " << std::endl;

										// ----------------------------------------------------------------
										// CONSTRUCTING SPECTRAL DENSITY MATRIX
										// ----------------------------------------------------------------

										if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
										{
											SpecDens *= 0.0;

											if (!ConstructSpecDensGeneral(1, ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
										}
										else
										{
											SpecDens *= 0.0;

											if (!ConstructSpecDensGeneral(0, ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
										}

#pragma omp parallel for num_threads(threads)
										for (l = 0; l < threads; l++)
										{
											*ptr_R[l] *= 0.0;
										}

#pragma omp parallel for firstprivate(tmp_R) shared(ampl_list, tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
										for (k = 0; k < num_op; k++)
										{
											// -----------------------------------------------------------------
											// CONSTRUCTING R MATRIX
											// -----------------------------------------------------------------

											tmp_R *= 0.0;

											if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[k]), SpecDens, tmp_R))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											*ptr_R[omp_get_thread_num()] += tmp_R;
										}
									}
								}
								else
								{
									if ((*interaction)->Properties()->Get("def_g", def_g) && def_g == 1)
									{
										this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;

#pragma omp parallel for num_threads(threads)
										for (l = 0; l < threads; l++)
										{
											*ptr_R[l] *= 0.0;
										}

// loop through all operators and construct R tensor in Liouville space
#pragma omp parallel for collapse(2) firstprivate(tmp_R, SpecDens, ampl_combined) shared(ampl_list, tau_c_list, domega, num_op, ptr_Tensors, ptr_R, interaction, def_specdens) num_threads(threads)
										for (k = 0; k < num_op; k++)
										{
											for (s = 0; s < num_op; s++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING SPECTRAL DENSITY MATRIX
												// ----------------------------------------------------------------

												ampl_combined = ampl_list[k] * ampl_list[s];

												if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
												{
													SpecDens *= 0.0;
													if (!ConstructSpecDensSpecific(1, ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
												}
												else
												{
													SpecDens *= 0.0;

													if (!ConstructSpecDensSpecific(0, ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
												}

												// ----------------------------------------------------------------
												// CONSTRUCTING R MATRIX
												// ----------------------------------------------------------------

												tmp_R *= 0.0;

												if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[s]), SpecDens, tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;

												if (k != s)
												{
													tmp_R *= 0.0;

													if (!Redfieldtensor((*ptr_Tensors[s]), (*ptr_Tensors[k]), SpecDens, tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;
												}
											}
										}
									}
									else
									{
										// loop through all operators and construct R tensor in Liouville space
										// ----------------------------------------------------------------
										// CONSTRUCTING SPECTRAL DENSITY MATRIX
										// ----------------------------------------------------------------

										this->Log() << "J is generally constructed for all operators - def_g == 0 " << std::endl;
										if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
										{
											SpecDens *= 0.0;

											if (!ConstructSpecDensGeneral(1, ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
										}
										else
										{
											SpecDens *= 0.0;

											if (!ConstructSpecDensGeneral(0, ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
										}

#pragma omp parallel for num_threads(threads)
										for (l = 0; l < threads; l++)
										{
											*ptr_R[l] *= 0.0;
										}

#pragma omp parallel for collapse(2) firstprivate(tmp_R) shared(ampl_list, tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
										for (k = 0; k < num_op; k++)
										{
											for (s = 0; s < num_op; s++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING R MATRIX
												// ----------------------------------------------------------------

												tmp_R *= 0.0;

												if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[s]), SpecDens, tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;

												if (k != s)
												{
													tmp_R *= 0.0;

													if (!Redfieldtensor((*ptr_Tensors[s]), (*ptr_Tensors[k]), SpecDens, tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}
													*ptr_R[omp_get_thread_num()] += tmp_R;
												}
											}
										}
									}
								}

								for (l = 0; l < threads; l++)
								{
									R += *ptr_R[l];
								}

								this->Log() << "Added relaxation matrix term for interaction " << (*interaction)->Name() << " of spin " << (*s1)->Name() << "." << std::endl;
							}
						}
						else
						{
							// Get groups with respect to interaction
							auto group1 = (*interaction)->Group1();
							auto group2 = (*interaction)->Group2();

							// Loop through groups to get all interaction [group1 = electrons, group2 = particles electrons interact with, commonly nuclei]
							for (auto s1 = group1.cbegin(); s1 < group1.cend(); s1++)
							{
								for (auto s2 = group2.cbegin(); s2 < group2.cend(); s2++)
								{
									if ((*interaction)->Properties()->Get("ops", ops) && ops == 1)
									{

										// --------------------------------------------------------
										// CREATION OF SPIN OPERATORS (Sz, Sx, Sy)
										// --------------------------------------------------------

										this->Log() << "Sx, Sy and Sz operator basis was chosen - ops == 1" << std::endl;

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sz()), *s1, *Sz1))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sx()), *s1, *Sx1))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s1)->Sy()), *s1, *Sy1))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s2)->Sz()), *s2, *Sz2))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s2)->Sx()), *s2, *Sx2))
										{
											return false;
										}

										if (!space.CreateOperator(arma::conv_to<arma::cx_mat>::from((*s2)->Sy()), *s2, *Sy2))
										{
											return false;
										}

										// Build double-spin operators
										*Sx1Sx2 = (*Sx1) * (*Sx2);
										*Sx1Sy2 = (*Sx1) * (*Sy2);
										*Sx1Sz2 = (*Sx1) * (*Sz2);
										*Sy1Sx2 = (*Sy1) * (*Sx2);
										*Sy1Sy2 = (*Sy1) * (*Sy2);
										*Sy1Sz2 = (*Sy1) * (*Sz2);
										*Sz1Sx2 = (*Sz1) * (*Sx2);
										*Sz1Sy2 = (*Sz1) * (*Sy2);
										*Sz1Sz2 = (*Sz1) * (*Sz2);

										// Put all tensors on pointer array
										num_op = 9;
										delete[] ptr_Tensors;
										ptr_Tensors = new arma::cx_mat *[num_op];

										ptr_Tensors[0] = Sx1Sx2;
										ptr_Tensors[1] = Sx1Sy2;
										ptr_Tensors[2] = Sx1Sz2;
										ptr_Tensors[3] = Sy1Sx2;
										ptr_Tensors[4] = Sy1Sy2;
										ptr_Tensors[5] = Sy1Sz2;
										ptr_Tensors[6] = Sz1Sx2;
										ptr_Tensors[7] = Sz1Sy2;
										ptr_Tensors[8] = Sz1Sz2;
									}
									else
									{
										// -------------------------------------------------------------------
										// Relaxation Matrix Construction for Interactions of Type::DoubleSpin
										// -------------------------------------------------------------------

										this->Log() << "Irreducible spherical tensor operator basis (Rank 0 & Rank 2) was chosen - ops == 0" << std::endl;

										// ----------------------------------------------------------------
										// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
										// ----------------------------------------------------------------

										if (!space.BlRk0TensorT0(*s1, *s2, *T0_rank_0))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=0
										if (!space.BlRk2SphericalTensorT0(*s1, *s2, *T0_rank_2))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=1
										if (!space.BlRk2SphericalTensorTp1(*s1, *s2, *Tp1))
										{
											this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=-1
										if (!space.BlRk2SphericalTensorTm1(*s1, *s2, *Tm1))
										{
											this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=2
										if (!space.BlRk2SphericalTensorTp2(*s1, *s2, *Tp2))
										{
											this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=-2
										if (!space.BlRk2SphericalTensorTm2(*s1, *s2, *Tm2))
										{
											this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
										num_op = 6;
										delete[] ptr_Tensors;
										ptr_Tensors = new arma::cx_mat *[num_op];
										ptr_Tensors[0] = T0_rank_0;
										ptr_Tensors[1] = T0_rank_2;
										ptr_Tensors[2] = Tm1;
										ptr_Tensors[3] = Tp1;
										ptr_Tensors[4] = Tm2;
										ptr_Tensors[5] = Tp2;

										// Construct spatial spherical Tensors
										SpinAPI::Tensor inTensor(0);
										if ((*interaction)->Properties()->Get("tensor", inTensor) && (*interaction)->Properties()->Get("coeff", coeff) && coeff == 1)
										{
											this->Log() << "Producing spatial shperical tensors with coupling tensor..." << std::endl;
											auto ATensor = (*interaction)->CouplingTensor();
											arma::cx_mat A(3, 3);
											arma::cx_vec Am(6);

											A.zeros();
											A = arma::conv_to<arma::cx_mat>::from((ATensor->LabFrame()));

											Am(0) = (1.0 / sqrt(3.0)) * (A(0, 0) + A(1, 1) + A(2, 2));
											Am(1) = (1.0 / sqrt(6.0)) * (3.0 * A(2, 2) - (A(0, 0) + A(1, 1) + A(2, 2)));
											Am(2) = 0.5 * (A(0, 2) + A(2, 0) - ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
											Am(3) = -0.5 * (A(0, 2) + A(2, 0) + ((arma::cx_double(0.0, 1.0)) * (A(1, 2) + A(2, 1))));
											Am(4) = 0.5 * (A(0, 0) - A(1, 1) - ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));
											Am(5) = 0.5 * (A(0, 0) - A(1, 1) + ((arma::cx_double(0.0, 1.0)) * (A(0, 1) + A(1, 0))));

											// Rank 0
											*ptr_Tensors[0] = Am(0) * (*ptr_Tensors[0]);
											// Rank 2
											*ptr_Tensors[1] = Am(1) * (*ptr_Tensors[1]);
											*ptr_Tensors[2] = -Am(3) * (*ptr_Tensors[2]);
											*ptr_Tensors[3] = -Am(2) * (*ptr_Tensors[3]);
											*ptr_Tensors[4] = Am(4) * (*ptr_Tensors[4]);
											*ptr_Tensors[5] = Am(5) * (*ptr_Tensors[5]);
										}
										else
										{
											// Norm
											*ptr_Tensors[0] = (*ptr_Tensors[0]);
											*ptr_Tensors[1] = (*ptr_Tensors[1]);
											*ptr_Tensors[2] = -(*ptr_Tensors[2]);
											*ptr_Tensors[3] = -(*ptr_Tensors[3]);
											*ptr_Tensors[4] = (*ptr_Tensors[4]);
											*ptr_Tensors[5] = (*ptr_Tensors[5]);
										}
									}

// Rotate tensors in eigenbasis of H0
#pragma omp parallel for firstprivate(ptr_Tensors) shared(eigen_vec) num_threads(threads)
									for (k = 0; k < num_op; k++)
									{
										*ptr_Tensors[k] = (eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec);
									}

									if ((*interaction)->Properties()->Get("slip", slip) && slip == 1)
									{
										this->Log() << "Calculating slippage of inital condition for interaction" << std::endl;
										arma::cx_mat rho0_new;
										rho0_new.set_size(size(rho0));
										rho0_new.zeros();

										if (!Slippage(ptr_Tensors, num_op, eig_val_mat, domega, rho0, static_cast<std::complex<double>>(tau_c_list[0]), rho0_new))
										{
											this->Log() << "There are problems with the construction of the intial condition for slippage - Please check your input." << std::endl;
											continue;
										}

										rho0 = rho0_new;
										rho0_new *= 0.0;
									}

									this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << std::endl;

									if ((*interaction)->Properties()->Get("terms", terms) && terms == 1)
									{
										this->Log() << "No cross-relaxation terms are requested" << std::endl;

										if ((*interaction)->Properties()->Get("def_g", def_g) && def_g == 1)
										{
											this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;

// loop through all operators and construct R tensor in Liouville space
#pragma omp parallel for num_threads(threads)
											for (int l = 0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

#pragma omp parallel for firstprivate(tmp_R, SpecDens, ampl_combined) shared(ampl_list, tau_c_list, domega, num_op, ptr_Tensors, ptr_R, interaction, def_specdens) num_threads(threads)
											for (k = 0; k < num_op; k++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING SPECTRAL DENSITY MATRIX
												// ----------------------------------------------------------------
												ampl_combined = ampl_list[k] * ampl_list[k];

												if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
												{
													SpecDens *= 0.0;

													if (!ConstructSpecDensSpecific(1, static_cast<std::complex<double>>(ampl_combined), static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
												}
												else
												{
													SpecDens *= 0.0;

													if (!ConstructSpecDensSpecific(0, ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
												}

												// -----------------------------------------------------------------
												// CONSTRUCTING R MATRIX
												// -----------------------------------------------------------------

												tmp_R *= 0.0;

												if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[k]), SpecDens, tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;
											}
										}
										else
										{
											this->Log() << "J is generally constructed for all operators - def_g == 0 " << std::endl;

											// ----------------------------------------------------------------
											// CONSTRUCTING SPECTRAL DENSITY MATRIX
											// ----------------------------------------------------------------

											if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
											{
												SpecDens *= 0.0;

												if (!ConstructSpecDensGeneral(1, ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
											}
											else
											{
												SpecDens *= 0.0;

												if (!ConstructSpecDensGeneral(0, ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
											}

#pragma omp parallel for num_threads(threads)
											for (l = 0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

#pragma omp parallel for firstprivate(tmp_R) shared(ampl_list, tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
											for (k = 0; k < num_op; k++)
											{
												// -----------------------------------------------------------------
												// CONSTRUCTING R MATRIX
												// -----------------------------------------------------------------

												tmp_R *= 0.0;

												if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[k]), SpecDens, tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;
											}
										}
									}
									else
									{
										if ((*interaction)->Properties()->Get("def_g", def_g) && def_g == 1)
										{
											this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;

#pragma omp parallel for num_threads(threads)
											for (l = 0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

// loop through all operators and construct R tensor in Liouville space
#pragma omp parallel for collapse(2) firstprivate(tmp_R, SpecDens, ampl_combined) shared(ampl_list, tau_c_list, domega, num_op, ptr_Tensors, ptr_R, interaction, def_specdens) num_threads(threads)
											for (k = 0; k < num_op; k++)
											{
												for (s = 0; s < num_op; s++)
												{
													// ----------------------------------------------------------------
													// CONSTRUCTING SPECTRAL DENSITY MATRIX
													// ----------------------------------------------------------------

													ampl_combined = ampl_list[k] * ampl_list[s];

													if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
													{
														SpecDens *= 0.0;
														if (!ConstructSpecDensSpecific(1, ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}

														SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
													}
													else
													{
														SpecDens *= 0.0;

														if (!ConstructSpecDensSpecific(0, ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}

														SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
													}

													// ----------------------------------------------------------------
													// CONSTRUCTING R MATRIX
													// ----------------------------------------------------------------

													tmp_R *= 0.0;

													if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[s]), SpecDens, tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;

														if (!Redfieldtensor((*ptr_Tensors[s]), (*ptr_Tensors[k]), SpecDens, tmp_R))
														{
															this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
															continue;
														}

														*ptr_R[omp_get_thread_num()] += tmp_R;
													}
												}
											}
										}
										else
										{
											// loop through all operators and construct R tensor in Liouville space
											// ----------------------------------------------------------------
											// CONSTRUCTING SPECTRAL DENSITY MATRIX
											// ----------------------------------------------------------------

											this->Log() << "J is generally constructed for all operators - def_g == 0 " << std::endl;
											if ((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
											{
												SpecDens *= 0.0;

												if (!ConstructSpecDensGeneral(1, ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
											}
											else
											{
												SpecDens *= 0.0;

												if (!ConstructSpecDensGeneral(0, ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= ((*interaction)->Prefactor() * (*interaction)->Prefactor());
											}

#pragma omp parallel for num_threads(threads)
											for (l = 0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

#pragma omp parallel for collapse(2) firstprivate(tmp_R) shared(ampl_list, tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
											for (k = 0; k < num_op; k++)
											{
												for (s = 0; s < num_op; s++)
												{
													// ----------------------------------------------------------------
													// CONSTRUCTING R MATRIX
													// ----------------------------------------------------------------

													tmp_R *= 0.0;

													if (!Redfieldtensor((*ptr_Tensors[k]), (*ptr_Tensors[s]), SpecDens, tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;

														if (!Redfieldtensor((*ptr_Tensors[s]), (*ptr_Tensors[k]), SpecDens, tmp_R))
														{
															this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
															continue;
														}
														*ptr_R[omp_get_thread_num()] += tmp_R;
													}
												}
											}
										}
									}

									for (l = 0; l < threads; l++)
									{
										R += *ptr_R[l];
									}

									this->Log() << "Added relaxation matrix term for interaction " << (*interaction)->Name() << " of spin " << (*s1)->Name() << "." << std::endl;
								}
							}
						}
					}
				}
				else
				{
					this->Log() << "No relaxation matrix term requested for interaction " << (*interaction)->Name() << "." << std::endl;
				}

				std::fill(ampl_list.begin(), ampl_list.end(), 0);
				std::fill(tau_c_list.begin(), tau_c_list.end(), 0);
				terms *= 0.0;
				SpecDens *= 0.0;
			}

			for (int l = 0; l < threads; l++)
			{
				delete ptr_R[l];
			}
			delete[] ptr_R;

			for (int l = 0; l < num_op; l++)
			{
				delete ptr_Tensors[l];
			}
			delete[] ptr_Tensors;

			// ---------------------------------------------------------------
			// SETUP COMPLETE HAMILTONIAN
			// ---------------------------------------------------------------
			// Transform  H0 into superspace
			arma::cx_mat lhs;
			arma::cx_mat rhs;
			arma::cx_mat H_SS;

			// Transforming into superspace
			space.SuperoperatorFromLeftOperator(eig_val_mat, lhs);
			space.SuperoperatorFromRightOperator(eig_val_mat, rhs);

			H_SS = lhs - rhs;

			// Get a matrix to collect all the terms (the total Liouvillian)
			arma::cx_mat A = arma::cx_double(0.0, -1.0) * H_SS;

			// Get the reaction operators, and add them to "A"
			arma::cx_mat K;

			if (!space.TotalReactionOperator(K))
			{
				this->Log() << "Warning: Failed to obtain matrix representation of the reaction operators!" << std::endl;
			}

			// Rotation into eigenbasis of H0
			K = eigen_vec.t() * K * eigen_vec;

			// Transform ReactionOperator into superspace
			arma::cx_mat Klhs;
			arma::cx_mat Krhs;
			arma::cx_mat K_SS;

			space.SuperoperatorFromLeftOperator(K, Klhs);
			space.SuperoperatorFromRightOperator(K, Krhs);

			K_SS = Klhs + Krhs;

			A -= K_SS;

			// Get the relaxation terms of other relaxation operators, assuming that they can just be added
			arma::cx_mat O_SS;

			for (auto t = (*i)->operators_cbegin(); t != (*i)->operators_cend(); t++)
			{
				space.UseSuperoperatorSpace(true);

				if (space.RelaxationOperatorFrameChange((*t), eigen_vec, O_SS))
				{

					A += O_SS;

					space.UseSuperoperatorSpace(false);
					this->Log() << "Added other relaxation operator \"" << (*t)->Name() << "\" to the Liouvillian.\n";
				}
				else
				{
					this->Log() << "There is a problem with operator \"" << (*t)->Name() << ". Please check.\n";
					space.UseSuperoperatorSpace(false);
				}
			}

			int secular = 0;
			this->Properties()->Get("secular", secular);

			if (secular == 1)
			{
				this->Log() << "Secular approximation is USED." << std::endl;

				double cutoff = 0.0;
				this->Properties()->Get("cutoff", cutoff);
				this->Log() << "Cutoff set to: " << cutoff << std::endl;

				int nn = domega.n_rows;

				std::cout << domega << std::endl;

				std::cout << "Elements of domega: " << nn << std::endl;
				std::cout << "Elements of R.n_col: " << R.n_cols << std::endl;
				std::cout << "Elements of R.n_row: " << R.n_rows << std::endl;

				int ra = 0;
				int rb = 0;
				int rc = 0;
				int rd = 0;
				int ii = 0;
				int jj = 0;
				int counter = 0;

				for (ra = 0; ra < nn; ra++)
				{
					for (rb = 0; rb < nn; rb++)
					{
						ii = ra * nn + rb;

						for (rc = 0; rc < nn; rc++)
						{
							for (rd = 0; rd < nn; rd++)
							{
								jj = rc * nn + rd;

								if (std::abs((domega(ra, rb) - domega(rc, rd))) >= cutoff)
								{
									R(ii, jj) *= 0.0;
									counter += 1;
								}
								else
								{
									std::cout << "domega(ra,rb): " << domega(ra, rb) << std::endl;
									std::cout << "domega(rc,rd): " << domega(rc, rd) << std::endl;
								}
							}
						}
					}
				}
				std::cout << "counter: " << counter << std::endl;
			}

			// Adding R tensor to whole hamiltonian
			A += R;

			// Transform density operator into superspace
			arma::cx_vec rho0vec;
			rho0vec *= 0.0;

			if (!space.OperatorToSuperspace(rho0, rho0vec))
			{
				this->Log() << "Failed to convert initial state density operator to superspace." << std::endl;
				continue;
			}

			// ---------------------------------------------------------------
			// DO PROPAGATION OF DENSITY OPERATOR
			// ---------------------------------------------------------------
			// Perform the calculation
			this->Log() << "Ready to perform calculation." << std::endl;
			arma::cx_vec result = solve(arma::conv_to<arma::cx_mat>::from(A), rho0vec);
			this->Log() << "Done with calculation." << std::endl;

			// Convert the resulting density operator back to its Hilbert space representation
			if (!space.OperatorFromSuperspace(result, rho0))
			{
				this->Log() << "Failed to convert resulting superspace-vector back to native Hilbert space." << std::endl;
				continue;
			}

			// ---------------------------------------------------------------
			// CALCULATE QUANTUM YIELDS
			// ---------------------------------------------------------------
			// Obtain the results

#pragma omp single
			{
				arma::cx_mat P;
				double sum_yield;
				sum_yield = 0.0;

				this->Data() << this->RunSettings()->CurrentStep() << " ";
				this->WriteStandardOutput(this->Data());

				// There are two result modes - either write results per transition or for each defined state
				if (this->productYieldsOnly)
				{
					// Loop through all defind transitions
					auto transitions = (*i)->Transitions();
					for (auto j = transitions.cbegin(); j < transitions.cend(); j++)
					{
						// Make sure that there is a state object
						if ((*j)->SourceState() == nullptr)
							continue;

						if (!space.GetState((*j)->SourceState(), P))
						{
							this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
							continue;
						}

						// Rotate projection operator in eigenbasis of H0
						P = (eigen_vec.t() * P * eigen_vec);

						// Return the yield for this transition
						this->Data() << (*j)->Rate() * std::abs(arma::trace(P * rho0)) << " ";
						sum_yield += ((*j)->Rate() * std::abs(arma::trace(P * rho0)));
					}
					this->Data() << sum_yield << " ";
					sum_yield *= 0.0;
					P *= 0.0;
				}
				else
				{
					// Loop through all states
					auto states = (*i)->States();
					for (auto j = states.cbegin(); j < states.cend(); j++)
					{
						if (!space.GetState((*j), P))
						{
							this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
							continue;
						}

						// Rotate projection operator in eigenbasis of H0
						P = (eigen_vec.t() * P * eigen_vec);

						// Return the yield for this state - note that no reaction rates are included here.
						this->Data() << std::abs(arma::trace(P * rho0)) << " ";
						sum_yield += (std::abs(arma::trace(P * rho0)));
					}

					this->Data() << sum_yield << " ";
					sum_yield *= 0.0;
					P *= 0.0;
				}

				this->Log() << "\nDone with SpinSystem \"" << (*i)->Name() << "\"" << std::endl;

				// Clean every variable just to be sure
				R *= 0.0;
				A *= 0.0;
				H *= 0.0;
				P *= 0.0;
				domega *= 0.0;
				H_SS *= 0.0;
				K_SS *= 0.0;
				eig_val_mat *= 0.0;
				eigen_val *= 0.0;
				eigen_vec *= 0.0;
				SpecDens *= 0.0;
			}
		}

		// Terminate the line in the data file after iteration through all spin systems
		this->Data() << std::endl;

		return true;
	}

	// -----------------------FUNCTIONS-----------------------------------------------

	// Writes the header of the data file (but can also be passed to other streams)
	void TaskStaticSSRedfield::WriteHeader(std::ostream &_stream)
	{
		_stream << "Step ";
		this->WriteStandardOutputHeader(_stream);

		// Get header for each spin system
		auto systems = this->SpinSystems();
		for (auto i = systems.cbegin(); i < systems.cend(); i++)
		{
			// Should yields be written per transition or per defined state?
			if (this->productYieldsOnly)
			{
				// Write each transition name
				auto transitions = (*i)->Transitions();
				for (auto j = transitions.cbegin(); j < transitions.cend(); j++)
					_stream << (*i)->Name() << "." << (*j)->Name() << ".yield ";
			}
			else
			{
				// Write each state name
				auto states = (*i)->States();
				for (auto j = states.cbegin(); j < states.cend(); j++)
					_stream << (*i)->Name() << "." << (*j)->Name() << " ";
			}
		}
		_stream << std::endl;
	}

	// Validation
	bool TaskStaticSSRedfield::Validate()
	{
		this->Properties()->Get("transitionyields", this->productYieldsOnly);

		return true;
	}

	// Construction Refield tensor
	bool TaskStaticSSRedfield::Redfieldtensor(const arma::cx_mat &_op1, const arma::cx_mat &_op2, const arma::cx_mat &_specdens, arma::cx_mat &_redfieldtensor)
	{

		arma::cx_mat one = arma::eye<arma::cx_mat>(arma::size(_specdens));

		_redfieldtensor *= 0.0;

		_redfieldtensor += arma::kron((_op1).t(), arma::conj((_op2).t() % (_specdens).st()));
		_redfieldtensor += arma::kron(((_op2).t() % (_specdens).st()), (_op1).st());
		_redfieldtensor -= arma::kron(_op1 * ((_op2).t() % (_specdens).st()), one);
		_redfieldtensor -= arma::kron(one, arma::conj(_op1 * (_op2.t() % _specdens.st())));

		// old - version

		// r = A1[a,c] * A2[d,b] * (S[c,a] + S[b,d])
		// _redfieldtensor += arma::kron(((_op1)%_specdens.st()),(_op2).st());
		// _redfieldtensor += arma::kron((_op1),(((_op2).st())%_specdens)); // (A2 * S.T).T = A2.T * S
		// if b == d: r -= sum(A1[a,:] * A2[:,c] * S[c,:])
		// _redfieldtensor -= arma::kron((_op1)*((_op2)%_specdens.st()),one);
		// if a == c: r -= sum(A1[d,:] * S[:,d] * A2[:,b])
		// _redfieldtensor -= arma::kron(one,(((_op1)%_specdens.st())*(_op2)).st());

		return true;
	}

	bool TaskStaticSSRedfield::ConstructSpecDensGeneral(const int &_spectral_function, const std::vector<double> &_ampl_list, const std::vector<double> &_tau_c_list, const arma::cx_mat &_domega, arma::cx_mat &_specdens)
	{
		if (_spectral_function == 1)
		{
// Solution  of spectral density : S = Ampl*(tau_c/(1+domega²*tauc²))
#pragma omp for
			for (auto ii = 0; ii < (int)_tau_c_list.size(); ii++)
			{
				_specdens += (static_cast<std::complex<double>>(_ampl_list[ii]) * (static_cast<std::complex<double>>(_tau_c_list[ii]) / (arma::cx_double(1.00, 0.00) + (pow(_domega, 2) * (pow(static_cast<std::complex<double>>(_tau_c_list[ii]), 2))))));
			}
		}

		if (_spectral_function == 0)
		{
// Solution of spectral density: S = Ampl/(1/tau_c - i * domega)
#pragma omp for
			for (auto ii = 0; ii < (int)_tau_c_list.size(); ii++)
			{
				_specdens += (static_cast<std::complex<double>>(_ampl_list[ii])) / ((1.00 / static_cast<std::complex<double>>(_tau_c_list[ii]) - (arma::cx_double(0.0, 1.0) * _domega)));
			}
		}

		return true;
	}

	bool TaskStaticSSRedfield::ConstructSpecDensSpecific(const int &_spectral_function, const std::complex<double> &_ampl, const std::complex<double> &_tau_c, const arma::cx_mat &_domega, arma::cx_mat &_specdens)
	{
		if (_spectral_function == 1)
		{
			// Solution  of spectral density : S = Ampl*(tau_c/(1+domega²*tauc²))
			_specdens = _ampl * (_tau_c / (arma::cx_double(1.00, 0.00) + (pow(_domega, 2) * (pow(_tau_c, 2)))));
		}

		if (_spectral_function == 0)
		{
			// Solution of spectral density: S = Ampl/(1/tau_c - i * domega)
			_specdens = _ampl / ((arma::cx_double(1.00, 0.00) / _tau_c) - (arma::cx_double(0.0, 1.0) * _domega));
		}

		return true;
	}

	bool TaskStaticSSRedfield::Slippage(arma::cx_mat **_ptr_Tensors, const int &_num_op, const arma::cx_mat &_eig_val_mat, const arma::cx_mat &_domega, const arma::cx_mat &_rho0, const std::complex<double> &_tau_c, arma::cx_mat &_rho0_new)
	{
		// Condition of Space: Hilbert
		// Input:
		// Pointer array of Tensors already in eigenbasis of HS: _ptr_Tensors
		// Number of Tensors for R tensor: _num_op
		// Eigenvalue diagonal matrix for HS: _eig_val_mat
		// Transition matrix omega for J's: _domega
		// Inital density matrix: rho0
		// Correlation time: _tau_c
		// Output:
		// Modified density matrix: rho0_new

		// Construct J with rank k
		arma::cx_mat Jk1, Jk2;

		// J^(k=1) = 1/ (i * domega + 1/tau_c)^(2)
		Jk1 = 1.0 / ((arma::cx_double(0.0, 1.0) * _domega + (1.0 / _tau_c)) * (arma::cx_double(0.0, 1.0) * _domega + (1.0 / _tau_c)));
		// J^(k=2) = 1/ (i * domega + 1/tau_c)^(3)
		Jk2 = 1.0 / ((arma::cx_double(0.0, 1.0) * _domega + (1.0 / _tau_c)) * (arma::cx_double(0.0, 1.0) * _domega + (1.0 / _tau_c)) * (arma::cx_double(0.0, 1.0) * _domega + (1.0 / _tau_c)));

		// Temporary matrix for summation
		arma::cx_mat tmp;
		tmp.set_size(size(Jk1));
		tmp.zeros();

		// Pointer arrays to store T_alpha operators
		arma::cx_mat *Tk1 = NULL;
		Tk1 = new arma::cx_mat[_num_op];
		arma::cx_mat *Tk2 = NULL;
		Tk2 = new arma::cx_mat[_num_op];

		// Contruct T_alpha operators

		for (int alpha = 0; alpha < _num_op; alpha++)
		{
			tmp *= 0.0;
			for (int beta = 0; beta < _num_op; beta++)
			{
				tmp += (Jk1 * (*_ptr_Tensors[beta]));
			}

			Tk1[alpha] = tmp;
		}

		for (int alpha = 0; alpha < _num_op; alpha++)
		{
			tmp *= 0.0;
			for (int beta = 0; beta < _num_op; beta++)
			{
				tmp += (Jk2 * (*_ptr_Tensors[beta]));
			}

			Tk2[alpha] = tmp;
		}

		// Construct Slippage operator
		arma::cx_mat slip_op;
		slip_op.set_size(size(Jk1));
		slip_op.zeros();

		for (int alpha = 0; alpha < _num_op; alpha++)
		{
			// T_alpha(k=1)
			slip_op += (*_ptr_Tensors[alpha]) * (Tk1[alpha]) * _rho0 - (Tk1[alpha]) * _rho0 * (*_ptr_Tensors[alpha]);

			// T_alpha(k=2)
			slip_op += (*_ptr_Tensors[alpha]) * (Tk2[alpha]) * _rho0 * (arma::cx_double(0.0, 1.0)) * _eig_val_mat - (*_ptr_Tensors[alpha]) * (Tk2[alpha]) * (arma::cx_double(0.0, 1.0)) * _eig_val_mat * _rho0;
			slip_op += -(Tk2[alpha]) * _rho0 * (arma::cx_double(0.0, 1.0)) * _eig_val_mat * (*_ptr_Tensors[alpha]) + (Tk2[alpha]) * (arma::cx_double(0.0, 1.0)) * _eig_val_mat * _rho0 * (*_ptr_Tensors[alpha]);
			slip_op += -(*_ptr_Tensors[alpha]) * (Tk2[alpha]) * _rho0 * (arma::cx_double(0.0, 1.0)) * _eig_val_mat + (Tk2[alpha]) * _rho0 * (*_ptr_Tensors[alpha]) * (arma::cx_double(0.0, 1.0)) * _eig_val_mat;
			slip_op += (arma::cx_double(0.0, 1.0)) * _eig_val_mat * (*_ptr_Tensors[alpha]) * (Tk2[alpha]) * _rho0 - (arma::cx_double(0.0, 1.0)) * _eig_val_mat * (Tk2[alpha]) * _rho0 * (*_ptr_Tensors[alpha]);

			// Hermitian conjugated (h.c.) -> .t() hermitian conjugated matrix
			slip_op += (*_ptr_Tensors[alpha]).t() * (Tk1[alpha]).t() * _rho0.t() - (Tk1[alpha]).t() * _rho0.t() * (*_ptr_Tensors[alpha]).t();

			slip_op += (*_ptr_Tensors[alpha]).t() * Tk2[alpha].t() * _rho0.t() * (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() - (*_ptr_Tensors[alpha]).t() * Tk2[alpha].t() * (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() * _rho0.t();
			slip_op += -Tk2[alpha].t() * _rho0.t() * (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() * (*_ptr_Tensors[alpha]).t() + Tk2[alpha].t() * (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() * _rho0.t() * (*_ptr_Tensors[alpha]).t();
			slip_op += -(*_ptr_Tensors[alpha]).t() * Tk2[alpha].t() * _rho0.t() * (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() + Tk2[alpha].t() * _rho0.t() * (*_ptr_Tensors[alpha]).t() * (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t();
			slip_op += (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() * (*_ptr_Tensors[alpha]).t() * Tk2[alpha].t() * _rho0.t() - (arma::cx_double(0.0, -1.0)) * _eig_val_mat.t() * Tk2[alpha].t() * _rho0.t() * (*_ptr_Tensors[alpha]).t();
		}

		// Update initial density matrix
		_rho0_new = _rho0 - slip_op;

		// Delete pointer arrays

		delete[] Tk1;
		delete[] Tk2;

		return true;
	}
}
