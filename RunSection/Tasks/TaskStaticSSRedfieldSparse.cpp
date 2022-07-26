//////////////////////////////////////////////////////////////////////////////
// MolSpin - Redfield Theory Task - developed by Luca Gerhards
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <omp.h>
#include "TaskStaticSSRedfieldSparse.h"
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
	TaskStaticSSRedfieldSparse::TaskStaticSSRedfieldSparse(const MSDParser::ObjectParser& _parser, const RunSection& _runsection)	: BasicTask(_parser,_runsection), reactionOperators(SpinAPI::ReactionOperatorType::Haberkorn), productYieldsOnly(false)
	{
		
	}
	
	TaskStaticSSRedfieldSparse::~TaskStaticSSRedfieldSparse()
	{
		
	}
	// -----------------------------------------------------
	// TaskStaticSSRedfield protected methods
	// -----------------------------------------------------
	bool TaskStaticSSRedfieldSparse::RunLocal()
	{
		this->Log() << "Running method StaticSSRedfieldSparse." << std::endl;
		
		// If this is the first step, write first part of header to the data file
		if(this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}
		
		// ----------------------------------------------------------------
		// SETTING UP SPIN SYSTEM AND GETTING DENSITY OPERATOR
		// ----------------------------------------------------------------
				
		// Loop through all SpinSystems
		auto systems = this->SpinSystems();
		for(auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Make sure we have an initial state
			auto initial_states = (*i)->InitialState();
			if(initial_states.size() < 1)
			{
				this->Log() << "Skipping SpinSystem \"" << (*i)->Name() << "\" as no initial state was specified." << std::endl;
				continue;
			}
			
			this->Log() << "\nStarting with SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
			
			// Obtain a SpinSpace to describe the system
			SpinAPI::SpinSpace space(*(*i));
			arma::sp_cx_mat rho0;
			space.UseSuperoperatorSpace(false);
			space.SetReactionOperatorType(this->reactionOperators);

			// Get the initial state
			for(auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
			{
				arma::sp_cx_mat tmp_rho0;
				if(!space.GetState(*j, tmp_rho0))
				{
					this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
					continue;
				}
				if(j == initial_states.cbegin())
					rho0 = tmp_rho0;
				else
					rho0 += tmp_rho0;
			}
			
			rho0 /= arma::trace(rho0);	// The density operator should have a trace of 1

			// ----------------------------------------------------------------
			// Get the Hamiltonian
			// ----------------------------------------------------------------

			arma::cx_mat H;
			if(!space.Hamiltonian(H))
			{
			this->Log() << "Failed to obtain Hamiltonian in superspace." << std::endl;
			continue;
			}

			// ----------------------------------------------------------------
			// DIAGONALIZATION OF H0// We need all of these operators
			// ----------------------------------------------------------------

			arma::cx_mat eigen_vec;		// To hold eigenvectors
			arma::vec eigen_val;	    // To hold eigenvalues

			this->Log() << "Starting diagonalization..." << std::endl;
			arma::eig_sym(eigen_val,eigen_vec,H);
			this->Log() << "Diagonalization done! Eigenvalues: " << eigen_val.n_elem << ", eigenvectors: " << eigen_vec.n_cols << std::endl;


			// ----------------------------------------------------------------
			//CONSTRUCTING TRANSITION MATRIX "domega" OUT OF EIGENVALUES OF H0
			// ----------------------------------------------------------------

			arma::sp_cx_mat domega(size(H));

			for (int k = 0; k < (eigen_val.size()); k++)
			{
				for (int s = 0; s < (eigen_val.size()); s++)
				{
						domega(k,s) = eigen_val(k) - eigen_val(s);
				}
			}

			// ---------------------------------------------------------------
			// SETUP RELAXATION OPERATOR
			// ---------------------------------------------------------------
			
			// ---------------------------------------------------------------
			// SETUP RELAXATION OPERATOR
			// ---------------------------------------------------------------
			// Redfield tensor
			arma::sp_cx_mat R;
			R.set_size(size(kron(H,H)));
			R.zeros();
			std::cout << "Size of R:" << R.size() << std::endl;

			// Temporary Redfield tensor
			arma::sp_cx_mat tmp_R;
			tmp_R.set_size(size(kron(H,H)));
			tmp_R.zeros();
			
			// R Tensor array pointer for parallelization
			arma::sp_cx_mat ** ptr_R = NULL;
			int threads;

			#pragma omp parallel 
			{
				threads = omp_get_num_threads();
			}

			ptr_R = new arma::sp_cx_mat* [threads];

			#pragma omp for					
			for (int l = 0; l < threads; l++) 
			{
   				ptr_R[l] = new arma::sp_cx_mat(tmp_R);
				*ptr_R[l] *= 0.0;   
			}


			// Spectral density matrix
			arma::sp_cx_mat SpecDens;
			SpecDens.set_size(size(domega));
			SpecDens.zeros();

			// Unit matrix
			arma::sp_cx_mat one;
			one.set_size(size(domega));
			one.eye();


			std::vector<double> tau_c_list;
			std::vector<double> ampl_list;


			// Defining variables for parameter of user
			int terms = 0;
			int define_ampl = 0;
			int ops = 0;
			int define_specdens = 0;

			// Defining varibales for loops
			int num_op;
			arma::sp_cx_mat ** ptr_Tensors = NULL;
			int k;
			int s;	

			// ------------------------------------------------------------------
			// STARTING WITH RELAXATION MATRIX CONSTRUCTION - GOOD LUCK
			// ------------------------------------------------------------------

			this->Log() << "Starting with construction of relaxation matrix." << std::endl;
			for(auto interaction = (*i)->interactions_cbegin(); interaction != (*i)->interactions_cend(); interaction++)
			{
				if((*interaction)->Properties()->GetList("tau_c", tau_c_list))
				{ 
					if((*interaction)->Properties()->GetList("g", ampl_list))
					{ 
						if((*interaction)->Type() == SpinAPI::InteractionType::SingleSpin)
						{
							// ------------------------------------------------------------------
							// Relaxation Matrix Construction for Interactions of Type::SingleSpin
							// ------------------------------------------------------------------
							
							// Get groups with respect to interaction
							auto group1 = (*interaction)->Group1();
																	
							// Loop through groups to get all interaction type
							for(auto s1 = group1.cbegin(); s1 != group1.cend(); s1++)
							{
								// Which operator basis was chosen:						
								if((*interaction)->Properties()->Get("ops", ops) && ops == 1)
								{
									// --------------------------------------------------------
								    // CREATION OF SPIN OPERATORS (Sz, Sx, Sy)
								    // --------------------------------------------------------

									this->Log() << "Sz, Sx and Sy operator basis was chosen - ops == 1" << std::endl;

									//Spin-Operators
									arma::sp_cx_mat *Sz = new arma::sp_cx_mat;
									arma::sp_cx_mat *Sx= new arma::sp_cx_mat;
									arma::sp_cx_mat *Sy= new arma::sp_cx_mat;

									if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s1)->Sz() ), *s1, *Sz))			
									{
										return false;
									}	

									if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s1)->Sx() ), *s1, *Sx))			
									{
										return false;
									}	

									if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s1)->Sy() ), *s1, *Sy))			
									{
										return false;
									}	
							
									// Put all tensors on pointer array
									num_op = 3;
									ptr_Tensors = new arma::sp_cx_mat * [num_op];
									ptr_Tensors[0] = Sx;
									ptr_Tensors[1] = Sy;
									ptr_Tensors[2] = Sz;

									// Rotate tensors in eigenbasis of H0
									#pragma omp parallel for
									for(k=0; k < num_op; k++)
									{
										*ptr_Tensors[k] = eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec;
									}
								}
								else
								{
									// ----------------------------------------------------------------
									// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
									// ----------------------------------------------------------------		

									this->Log() << "Irreducible spherical tensor operator basis (Rank 0 & Rank 2) was chosen - ops == 0" << std::endl;
									// T0 for rank 0 & 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - space.Rk1SphericalTensorXXX)
									arma::sp_cx_mat *T0_rank_0 = new arma::sp_cx_mat;
									arma::sp_cx_mat *T0_rank_2 = new arma::sp_cx_mat;
									
									// Tp1 & T1m for rank 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - space.Rk1SphericalTensorXXX)
									arma::sp_cx_mat *Tp1 = new arma::sp_cx_mat;
									arma::sp_cx_mat *Tm1 = new arma::sp_cx_mat;

									// Tp2 & Tp2 for rank 2
									arma::sp_cx_mat *Tp2 = new arma::sp_cx_mat;
									arma::sp_cx_mat *Tm2 = new arma::sp_cx_mat;

									this->Log() << "Using magentic field to construct SingleSpin irreducible tensors." << std::endl;

									arma::vec static_field = (*interaction)->Field();
									arma::cx_vec complex_conv_static_field;
									complex_conv_static_field.set_size(size(static_field));
									
									complex_conv_static_field[0] = (static_field[0] + (arma::cx_double(0.0, 1.0) * static_field[1]));
									complex_conv_static_field[1] = (static_field[0] - (arma::cx_double(0.0, 1.0) * static_field[1]));
									complex_conv_static_field[2] = static_cast<std::complex<double>>(static_field[2]);
									
									// Rank 0 tensor with m=0
									*T0_rank_0 = -(1/sqrt(3.0))*one;

									// Rank 2 tensor with m=0
									if(!space.LRk2SphericalTensorT0(*s1, complex_conv_static_field, *T0_rank_2))
									{
										this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
										continue;
									}
																				
									// Rank 2 tensor with m=1
									if(!space.LRk2SphericalTensorTp1(*s1, complex_conv_static_field, *Tp1))
									{
										this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
										continue;
									} 

									// Rank 2 tensor with m=-1
									if(!space.LRk2SphericalTensorTm1(*s1, complex_conv_static_field, *Tm1))
									{
										this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
										continue;
									}
										
									// Rank 2 tensor with m=2
									if(!space.LRk2SphericalTensorTp2(*s1, complex_conv_static_field, *Tp2))
									{
										this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and Field" << "! Skipping." << std::endl;
										continue;
									}
																		
									// Rank 2 tensor with m=-2
									if(!space.LRk2SphericalTensorTm2(*s1, complex_conv_static_field, *Tm2))
									{
										this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and Field" << "! Skipping." << std::endl;
										continue;
									}
									// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
									num_op = 6;
									ptr_Tensors = new arma::sp_cx_mat* [num_op];							
									ptr_Tensors[0] = T0_rank_0;
									ptr_Tensors[1] = T0_rank_2; 
									ptr_Tensors[2] = Tm1; 
									ptr_Tensors[3] = Tp1;
									ptr_Tensors[4] = Tm2;
									ptr_Tensors[5] = Tp2;

									// Rotate tensors in eigenbasis of H0
									#pragma omp parallel for
									for(k=0; k < num_op; k++)
									{
										*ptr_Tensors[k] = eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec;
									}
								}

								this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << std::endl;

								if((*interaction)->Properties()->Get("terms", terms ) && terms == 1)
								{
									if((*interaction)->Properties()->Get("define_ampl", define_ampl ) && define_ampl == 1)
									{
										this->Log() << "No cross-relaxation terms are requested" << std::endl;
										this->Log() << "Setting J up for each operator separatly - define_ampl == 1 " << std::endl;

										// loop through all operators and construct R tensor in Liouville space
										#pragma omp parallel for firstprivate(tmp_R, SpecDens) shared(ampl_list,tau_c_list, domega)
										for(k=0; k < num_op;  k++)
										{
											// ----------------------------------------------------------------
											// CONSTRUCTING SPECTRAL DENSITY MATRIX
											// ----------------------------------------------------------------

											if((*interaction)->Properties()->Get("define_specdens", define_specdens) && define_specdens == 1)
											{
												SpecDens *= 0.0;
												if(!ConstructSpecDensSpecificSparse(1,static_cast<std::complex<double>> (ampl_list[k]), static_cast<std::complex<double>>(tau_c_list[k]), domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}
											else
											{
												SpecDens *= 0.0;
												
												if(!ConstructSpecDensSpecificSparse(0,static_cast<std::complex<double>> (ampl_list[k]), static_cast<std::complex<double>>(tau_c_list[k]), domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}

											// -----------------------------------------------------------------
											// CONSTRUCTING R MATRIX 
											// -----------------------------------------------------------------
											tmp_R *= 0.0;
											*ptr_R[omp_get_thread_num()] *= 0.0;

											if(!RedfieldtensorSparse((*ptr_Tensors[k]),(*ptr_Tensors[k]),SpecDens,tmp_R))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											*ptr_R[omp_get_thread_num()] += tmp_R;
										}
	
										#pragma omp for
										for (int l=0; l < threads; l++)
										{
											R += *ptr_R[l];
										}
									
									}
									else
									{
										this->Log() << "J is generally constructed for all operators - define_ampl == 0 " << std::endl;
										
										// ----------------------------------------------------------------
										// CONSTRUCTING SPECTRAL DENSITY MATRIX
										// ----------------------------------------------------------------	
										if((*interaction)->Properties()->Get("define_specdens", define_specdens) && define_specdens == 1)
										{
											SpecDens *= 0.0;

											if(!ConstructSpecDensGeneralSparse(1,ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= (*interaction)->Prefactor();
										}
										else
										{
											SpecDens *= 0.0;

											if(!ConstructSpecDensGeneralSparse(0,ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

											SpecDens *= (*interaction)->Prefactor();
										}
										
										#pragma omp parallel for firstprivate(tmp_R) shared(ampl_list,tau_c_list, domega, SpecDens)
										for(k=0; k < num_op;  k++)
										{						
											// -----------------------------------------------------------------
											// CONSTRUCTING R MATRIX 
											// -----------------------------------------------------------------
											tmp_R *= 0.0;
											*ptr_R[omp_get_thread_num()] *= 0.0;
										
											if(!RedfieldtensorSparse((*ptr_Tensors[k]),(*ptr_Tensors[k]),SpecDens,tmp_R))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											*ptr_R[omp_get_thread_num()] += tmp_R;
										}
										
										#pragma omp for
										for (int l=0; l < threads; l++)
										{
											R += *ptr_R[l];
										}
									}
								}	
								else
								{
									if((*interaction)->Properties()->Get("define_ampl", define_ampl ) && define_ampl == 1)
									{
										this->Log() << "Setting J up for each operator separatly - define_ampl == 1 " << std::endl;
										
										// loop through all operators and construct R tensor in Liouville space	
										#pragma omp parallel for collapse(2) firstprivate(tmp_R,SpecDens) shared(ampl_list,tau_c_list, domega)
										for(k=0; k < num_op;  k++)
										{
											for(s=0; s < num_op; s++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING SPECTRAL DENSITY MATRIX
												// ----------------------------------------------------------------
												if((*interaction)->Properties()->Get("define_specdens", define_specdens) && define_specdens == 1)
												{
													SpecDens *= 0.0;
													if(!ConstructSpecDensSpecificSparse(1,static_cast<std::complex<double>> (ampl_list[k]), static_cast<std::complex<double>>(tau_c_list[k]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= (*interaction)->Prefactor();
												}
												else
												{
													SpecDens *= 0.0;
													
													if(!ConstructSpecDensSpecificSparse(0,static_cast<std::complex<double>> (ampl_list[k]), static_cast<std::complex<double>>(tau_c_list[k]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}
													
													SpecDens *= (*interaction)->Prefactor();
												}
												
												// ----------------------------------------------------------------
												// CONSTRUCTING R MATRIX 
												// ----------------------------------------------------------------
												tmp_R *= 0.0;
												*ptr_R[omp_get_thread_num()] *= 0.0;

												if(!RedfieldtensorSparse((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;

												if (k != s)
												{
													tmp_R *= 0.0;

													if(!RedfieldtensorSparse((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}
												}										

												*ptr_R[omp_get_thread_num()] += tmp_R;
											}
										}			
										
										#pragma omp for
										for (int l=0; l < threads; l++)
										{
											R += *ptr_R[l];
										}
									}
									else
									{
										// loop through all operators and construct R tensor in Liouville space	
										// ----------------------------------------------------------------
										// CONSTRUCTING SPECTRAL DENSITY MATRIX
										// ----------------------------------------------------------------
										this->Log() << "J is generally constructed for all operators - define_ampl == 0 " << std::endl;
										if((*interaction)->Properties()->Get("define_specdens", define_specdens) && define_specdens == 1)
										{
											SpecDens *= 0.0;
										
											if(!ConstructSpecDensGeneralSparse(1,ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= (*interaction)->Prefactor();
										}
										else
										{
											SpecDens *= 0.0;

											if(!ConstructSpecDensGeneralSparse(0,ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= (*interaction)->Prefactor();
										}
										
										#pragma omp parallel for collapse(2) firstprivate(tmp_R) shared(ampl_list,tau_c_list, domega, SpecDens)
										for(k=0; k < num_op;  k++)
										{
											for(s=0; s < num_op; s++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING R MATRIX 
												// ----------------------------------------------------------------
												tmp_R *= 0.0;
												*ptr_R[omp_get_thread_num()] *= 0.0;
											
												if(!RedfieldtensorSparse((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;

												if (k != s)
												{
													tmp_R *= 0.0;

													if(!RedfieldtensorSparse((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}
												}
												
												*ptr_R[omp_get_thread_num()] += tmp_R;
											}

										}

										#pragma omp for
										for (int l=0; l < threads; l++)
										{
											R += *ptr_R[l];
										}
										
									}
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
							for(auto s1 = group1.cbegin(); s1 != group1.cend(); s1++)
							{
								for(auto s2 = group2.cbegin(); s2 != group2.cend(); s2++)
								{
									if((*interaction)->Properties()->Get("ops", ops) && ops == 1)
									{

										// --------------------------------------------------------
										// CREATION OF SPIN OPERATORS (Sz, Sx, Sy)
										// --------------------------------------------------------

										this->Log() << "Sz, Sx and Sy operator basis was chosen - ops == 1" << std::endl;

										//Spin-Operators
										arma::sp_cx_mat *Sz1 = new arma::sp_cx_mat;
										arma::sp_cx_mat *Sx1= new arma::sp_cx_mat;
										arma::sp_cx_mat *Sy1= new arma::sp_cx_mat;
										arma::sp_cx_mat *Sz2 = new arma::sp_cx_mat;
										arma::sp_cx_mat *Sx2= new arma::sp_cx_mat;
										arma::sp_cx_mat *Sy2= new arma::sp_cx_mat;

										if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s1)->Sz() ), *s1, *Sz1))			
										{
											return false;
										}	

										if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s1)->Sx() ), *s1, *Sx1))			
										{
											return false;
										}	

										if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s1)->Sy() ), *s1, *Sy1))			
										{
											return false;
										}	
										
										if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s2)->Sz() ), *s1, *Sz2))			
										{
											return false;
										}	

										if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s2)->Sx() ), *s1, *Sx2))			
										{
											return false;
										}	

										if(!space.CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( (*s2)->Sy() ), *s1, *Sy2))			
										{
											return false;
										}	

										// Put all tensors on pointer array
										num_op = 6;
										ptr_Tensors = new arma::sp_cx_mat * [num_op];
										ptr_Tensors[0] = Sx1;
										ptr_Tensors[1] = Sy1;
										ptr_Tensors[2] = Sz1;
										ptr_Tensors[3] = Sx2;
										ptr_Tensors[4] = Sy2;
										ptr_Tensors[5] = Sz2;

										// Rotate tensors in eigenbasis of H0
										#pragma omp parallel for
										for(k=0; k < num_op; k++)
										{
												*ptr_Tensors[k] = eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec;
										}
									
									}
									else
									{
										// -------------------------------------------------------------------
										// Relaxation Matrix Construction for Interactions of Type::DoubleSpin
										// -------------------------------------------------------------------

										this->Log() << "Irreducible spherical tensor operator basis (Rank 0 & Rank 2) was chosen - ops == 0" << std::endl;
										
										// T0 for rank 0 & 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - space.Rk1SphericalTensorXXX)
										arma::sp_cx_mat *T0_rank_0 = new arma::sp_cx_mat;
										arma::sp_cx_mat *T0_rank_2 = new arma::sp_cx_mat;
												
										// Tp1 & T1m for rank 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - space.Rk1SphericalTensorXXX)
										arma::sp_cx_mat *Tp1 = new arma::sp_cx_mat;
										arma::sp_cx_mat *Tm1 = new arma::sp_cx_mat;

										// Tp2 & Tp2 for rank 2
										arma::sp_cx_mat *Tp2 = new arma::sp_cx_mat;
										arma::sp_cx_mat *Tm2 = new arma::sp_cx_mat;
					
										// ----------------------------------------------------------------------
										// TODO: Sometimes the Coupling tensors are included in relaxation
										// ----------------------------------------------------------------------
										// auto ATensor = (*interaction)->CouplingTensor();
										// auto A = ATensor->LabFrame();
										// a0 = -(A(1,1)+A(2,2)+A(3,3));
										// Ax = (2.0*A(3,3)-((A(1,1)+A(2,2))));
										// Rh = (A(1,1)-A(2,2));
										// ----------------------------------------------------------------------

										// ----------------------------------------------------------------
										// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
										// ----------------------------------------------------------------								
										// Rank 0 tensor with m=0
										*T0_rank_0 = -(1/sqrt(3.0))*one;

										// Rank 2 tensor with m=0
										if(!space.BlRk2SphericalTensorT0(*s1, *s2, *T0_rank_2))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}
																					
										// Rank 2 tensor with m=1
										if(!space.BlRk2SphericalTensorTp1(*s1, *s2, *Tp1))
										{
											this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										} 

										// Rank 2 tensor with m=-1
										if(!space.BlRk2SphericalTensorTm1(*s1, *s2, *Tm1))
										{
											this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}
													
										// Rank 2 tensor with m=2
										if(!space.BlRk2SphericalTensorTp2(*s1, *s2, *Tp2))
										{
											this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}
																					
										// Rank 2 tensor with m=-2
										if(!space.BlRk2SphericalTensorTm2(*s1, *s2, *Tm2))
										{
											this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
											continue;
										}

										// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
										num_op = 6;						

										ptr_Tensors = new arma::sp_cx_mat* [num_op];
										ptr_Tensors[0] = T0_rank_0;
										ptr_Tensors[1] = T0_rank_2; 
										ptr_Tensors[2] = Tm1; 
										ptr_Tensors[3] = Tp1;
										ptr_Tensors[4] = Tm2;
										ptr_Tensors[5] = Tp2;

										// Rotate tensors in eigenbasis of H0
										#pragma omp parallel for
										for(k=0; k < num_op; k++)
										{
											*ptr_Tensors[k] = eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec;
										}
									}		

									this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << " and spin " << (*s2)->Name() << std::endl;

									if((*interaction)->Properties()->Get("define_ampl", define_ampl ) && define_ampl == 1)
									{
										this->Log() << "Setting J up for each operator separatly - define_ampl == 1 " << std::endl;

										#pragma omp parallel for collapse(2) firstprivate(tmp_R,SpecDens) shared(ampl_list,tau_c_list, domega)
										for(k=0; k < num_op;  k++)
										{
											for(s=0; s < num_op; s++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING SPECTRAL DENSITY MATRIX
												// ----------------------------------------------------------------
												if((*interaction)->Properties()->Get("define_specdens", define_specdens) && define_specdens == 1)
												{
													SpecDens *= 0.0;
													if(!ConstructSpecDensSpecificSparse(1,static_cast<std::complex<double>> (ampl_list[k]), static_cast<std::complex<double>>(tau_c_list[k]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= (*interaction)->Prefactor();

												}
												else
												{
													SpecDens *= 0.0;
															
													if(!ConstructSpecDensSpecificSparse(0,static_cast<std::complex<double>> (ampl_list[k]), static_cast<std::complex<double>>(tau_c_list[k]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}
															
													SpecDens *= (*interaction)->Prefactor();
												}

												// -----------------------------------------------------------------------------------------
												// CONSTRUCTING R MATRIX 
												// -----------------------------------------------------------------------------------------
												tmp_R *= 0.0;
												*ptr_R[omp_get_thread_num()] *= 0.0;

												if(!RedfieldtensorSparse((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;

												if (k != s)
												{
													tmp_R *= 0.0;
													if(!RedfieldtensorSparse((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;
											}	
										}
									}
									else
									{
										// ----------------------------------------------------------------
										// CONSTRUCTING SPECTRAL DENSITY MATRIX
										// ----------------------------------------------------------------

										if((*interaction)->Properties()->Get("define_specdens", define_specdens) && define_specdens == 1)
										{
											SpecDens *= 0.0;

											if(!ConstructSpecDensGeneralSparse(1,ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= (*interaction)->Prefactor();
										}
										else
										{
											SpecDens *= 0.0;

											if(!ConstructSpecDensGeneralSparse(0,ampl_list, tau_c_list, domega, SpecDens))
											{
												this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
												continue;
											}

											SpecDens *= (*interaction)->Prefactor();
										}

										#pragma omp parallel for collapse(2) firstprivate(tmp_R) shared(ampl_list,tau_c_list, domega, SpecDens)
										for(k=0; k < num_op;  k++)
										{
											for(s=0; s < num_op; s++)
											{
												// -----------------------------------------------------------------------------------------
												// CONSTRUCTING R MATRIX 
												// -----------------------------------------------------------------------------------------
												tmp_R *= 0.0;
												*ptr_R[omp_get_thread_num()] *= 0.0;

												if(!RedfieldtensorSparse((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;

												if (k != s)
												{
													tmp_R *= 0.0;
													if(!RedfieldtensorSparse((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}
												}

												*ptr_R[omp_get_thread_num()] += tmp_R;
											}	
										}
									}
									
									#pragma omp for
									for (int l=0; l < threads; l++)
									{
										R += *ptr_R[l];
									}

									this->Log() << "Added relaxation matrix term for interaction " << (*interaction)->Name() << " between spins " << (*s1)->Name() << " and " << (*s2)->Name() << "." << std::endl;
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
				SpecDens *=0.0;
				delete ptr_Tensors;
			}

			delete ptr_R;
			// ---------------------------------------------------------------

			// Add 1/(planck bar)²
			R *= 8.794e+1;

			// ---------------------------------------------------------------
			// SETUP COMPLETE HAMILTONIAN
			// ---------------------------------------------------------------

			// Transform  H0 into superspace
			arma::sp_cx_mat lhs;
			arma::sp_cx_mat rhs;
			arma::sp_cx_mat eig_val_mat;
			arma::sp_cx_mat H_SS;

			// Constructing diagonal matrix with eigenvalues of H0
			eig_val_mat = arma::diagmat(arma::conv_to<arma::cx_mat>::from(eigen_val));

			// Transforming into superspace
			space.SuperoperatorFromLeftOperator(eig_val_mat, lhs);
			space.SuperoperatorFromRightOperator(eig_val_mat, rhs);

		   	H_SS = lhs - rhs;

		        // Get a matrix to collect all the terms (the total Liouvillian)
			arma::sp_cx_mat A = arma::cx_double(0.0, -1.0) * H_SS;
			
			// Get the reaction operators, and add them to "A"
			arma::sp_cx_mat K;

			if(!space.TotalReactionOperator(K))
			{
				this->Log() << "Warning: Failed to obtain matrix representation of the reaction operators!" << std::endl;
			}

			// Rotation into eigenbasis of H0
			K = eigen_vec.t() * K * eigen_vec;

			// Transform ReactionOperator into superspace
			arma::sp_cx_mat Klhs;
			arma::sp_cx_mat Krhs;
			arma::sp_cx_mat K_SS;
			
			space.SuperoperatorFromLeftOperator(K, Klhs);
			space.SuperoperatorFromRightOperator(K, Krhs);

		        K_SS = Klhs + Krhs;

			A -= K_SS;

			// Adding R tensor to whole hamiltonian
			A += R;

			// Rotate density operator in eigenbasis of H0
			rho0 = (eigen_vec.t() * rho0 * eigen_vec);
						
			// Transform density operator into superspace
			arma::cx_vec rho0vec;
			rho0vec *= 0.0;

			if(!space.OperatorToSuperspace(rho0, rho0vec))
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
			if(!space.OperatorFromSuperspace(result, rho0))
			{
				this->Log() << "Failed to convert resulting superspace-vector back to native Hilbert space." << std::endl;
				continue;
			}

			// ---------------------------------------------------------------
			// CALCULATE QUANTUM YIELDS
			// ---------------------------------------------------------------

			// Obtain the results
			arma::cx_mat P;
			double sum_yield;
			this->Data() << this->RunSettings()->CurrentStep() << " ";
			this->WriteStandardOutput(this->Data());
			
			// There are two result modes - either write results per transition or for each defined state
			if(this->productYieldsOnly)
			{
				// Loop through all defind transitions
				auto transitions = (*i)->Transitions();
				for(auto j = transitions.cbegin(); j != transitions.cend(); j++)
				{
					// Make sure that there is a state object
					if((*j)->SourceState() == nullptr)
						continue;
					
					if(!space.GetState((*j)->SourceState(), P))
					{
						this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
						continue;
					}
					
					// Rotate projection operator in eigenbasis of H0
					P = (eigen_vec.t() * P * eigen_vec);

					sum_yield += (*j)->Rate() * std::abs(arma::trace(P * rho0));

					// Return the yield for this transition
					this->Data() << (*j)->Rate() * std::abs(arma::trace(P * rho0)) << " ";


				}
				this->Data()<< sum_yield << " ";
                                sum_yield *= 0.0;
                                P *= 0.0;
			}
			else
			{
				// Loop through all states
				auto states = (*i)->States();
				for(auto j = states.cbegin(); j != states.cend(); j++)
				{
					if(!space.GetState((*j), P))
					{
						this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
						continue;
					}
					
					// Rotate projection operator in eigenbasis of H0
					P = (eigen_vec.t() * P * eigen_vec);

					// Return the yield for this state - note that no reaction rates are included here.
					this->Data() << std::abs(arma::trace(P * rho0)) << " ";
				}
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
		
		// Terminate the line in the data file after iteration through all spin systems
		this->Data() << std::endl;
		
		return true;
	}

// --------------------------------------------------------------------------------------------------------------------------------------

	// Writes the header of the data file (but can also be passed to other streams)
	void TaskStaticSSRedfieldSparse::WriteHeader(std::ostream& _stream)
	{
		_stream << "Step ";
		this->WriteStandardOutputHeader(_stream);
		
		// Get header for each spin system
		auto systems = this->SpinSystems();
		for(auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Should yields be written per transition or per defined state?
			if(this->productYieldsOnly)
			{
				// Write each transition name
				auto transitions = (*i)->Transitions();
				for(auto j = transitions.cbegin(); j != transitions.cend(); j++)
					_stream << (*i)->Name() << "." << (*j)->Name() << ".yield ";
			}
			else
			{
				// Write each state name
				auto states = (*i)->States();
				for(auto j = states.cbegin(); j != states.cend(); j++)
					_stream << (*i)->Name() << "." << (*j)->Name() << " ";
			}
		}
		_stream << std::endl;
	}
	
	// Validation
	bool TaskStaticSSRedfieldSparse::Validate()
	{
		this->Properties()->Get("transitionyields",this->productYieldsOnly);
		
		return true;
	}

	//Construction Refield tensor
	bool TaskStaticSSRedfieldSparse::RedfieldtensorSparse(const arma::sp_cx_mat& _op1, const arma::sp_cx_mat& _op2, const arma::sp_cx_mat& _specdens, arma::sp_cx_mat& _redfieldtensor) 
	{

		arma::sp_cx_mat one = arma::eye<arma::sp_cx_mat>(arma::size(_specdens));

		_redfieldtensor*=0.0;

		// r = A1[a,c] * A2[d,b] * (S[c,a] + S[b,d])
		_redfieldtensor += arma::kron((_op1)%_specdens.st(),(_op2).st());
		_redfieldtensor += arma::kron((_op1),(((_op2).st())%_specdens)); // (A2 * S.T).T = A2.T * S
	    // if b == d: r -= sum(A1[a,:] * A2[:,c] * S[c,:])
		_redfieldtensor -= arma::kron((_op1)*((_op2)%_specdens.st()),one);
		// if a == c: r -= sum(A1[d,:] * S[:,d] * A2[:,b])
		_redfieldtensor -= arma::kron(one,(((_op1)%_specdens.st())*(_op2)).st());

		// --------------------------------------------------------------------------------------------
		// Algorithm provided by Dr. Daniel Kattnig, University of Exeter, UK
		// r = A1[a,c] * A2[d,b] * (S[c,a] + S[b,d])
		// R  = kron((*ptr_Tensors[k])%SpecDens.st(),(*ptr_Tensors[s]).st());
		// R += kron((*ptr_Tensors[k]),((*ptr_Tensors[s]).st()%SpecDens)); // (A2 * S.T).T = A2.T * S
		// if b == d: r -= sum(A1[a,:] * A2[:,c] * S[c,:])
		// R -= kron((*ptr_Tensors[k])*((*ptr_Tensors[s])%SpecDens.st()),one);
		// if a == c: r -= sum(A1[d,:] * S[:,d] * A2[:,b])
		// R -= kron(one,(((*ptr_Tensors[k])%SpecDens.st())*(*ptr_Tensors[s])).st());	

		//if(k!=s){same routine but switch indices k and s in equation} 
		// --------------------------------------------------------------------------------------------
		return true;

	}

	bool TaskStaticSSRedfieldSparse::ConstructSpecDensGeneralSparse(const int& _spectral_function, const std::vector<double>& _ampl_list, const std::vector<double>& _tau_c_list,const arma::sp_cx_mat& _domega, arma::sp_cx_mat& _specdens)
	{
		if(_spectral_function == 0)
		{
			// Solution  of spectral density : S = prefactor*Ampl*(tau_c/(1+domega²*tauc²))
			#pragma omp for
			for (auto ii = 0; ii != _tau_c_list.size(); ii++)
			{
			   	_specdens += ((static_cast<std::complex<double>> (_ampl_list[ii])) * static_cast<std::complex<double>>(_tau_c_list[ii]))/(1.00+(_domega*_domega) * (static_cast<std::complex<double>>(_tau_c_list[ii]) * static_cast<std::complex<double>>(_tau_c_list[ii]))); 
			}
		}

		if(_spectral_function == 1)
		{
			// Solution of spectral density: S = prefactor*Ampl/(1/tau_c - i * domega)
			#pragma omp for
			for (auto ii = 0; ii != _tau_c_list.size(); ii++)
				{
					_specdens +=  (static_cast<std::complex<double>>(_ampl_list[ii])) / ((1.00/static_cast<std::complex<double>>(_tau_c_list[ii]) - (arma::cx_double(0.0, 1.0)*_domega)));
				}	
		}

		return true;
	}

	bool TaskStaticSSRedfieldSparse::ConstructSpecDensSpecificSparse(const int& _spectral_function, const std::complex<double>& _ampl, const std::complex<double>& _tau_c,const arma::sp_cx_mat& _domega, arma::sp_cx_mat& _specdens)
	{
		if(_spectral_function == 0)
		{
			// Solution  of spectral density : S = prefactor*Ampl*(tau_c/(1+domega²*tauc²))
			_specdens = (_ampl * _tau_c)/(1.00+(_domega * _domega * _tau_c * _tau_c)); 

		}

		if(_spectral_function == 1)
		{
			// Solution of spectral density: S = prefactor*Ampl/(1/tau_c - i * domega)
			_specdens = _ampl / ((1.00/_tau_c) - (arma::cx_double(0.0, 1.0)*_domega));
		}

		return true;

	}

}

