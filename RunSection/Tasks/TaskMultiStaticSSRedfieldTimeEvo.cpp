/////////////////////////////////////////////////////////////////////////
// TaskMultiStaticSSRedfieldTimeEvo implementation (RunSection module)

// -- Multi-system version: Allows transitions between SpinSystems --
// 
// Molecular Spin Dynamics Software - developed by Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TaskMultiStaticSSRedfieldTimeEvo.h"
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "ObjectParser.h"

#include <iostream>
#include <omp.h>
#include <memory>
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "Spin.h"
#include "Interaction.h"
#include "Operator.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskMultiStaticSSRedfieldTimeEvo Constructors and Destructor
	// -----------------------------------------------------
	TaskMultiStaticSSRedfieldTimeEvo::TaskMultiStaticSSRedfieldTimeEvo(const MSDParser::ObjectParser& _parser, const RunSection& _runsection)	: BasicTask(_parser,_runsection), timestep(1.0), totaltime(1.0e+4),
																																	reactionOperators(SpinAPI::ReactionOperatorType::Haberkorn)
	{
		
	}
	
	TaskMultiStaticSSRedfieldTimeEvo::~TaskMultiStaticSSRedfieldTimeEvo()
	{
		
	}
	// -----------------------------------------------------
	// TaskMultiStaticSSRedfieldTimeEvo protected methods
	// -----------------------------------------------------
	bool TaskMultiStaticSSRedfieldTimeEvo::RunLocal()
	{
		this->Log() << "Running method StaticSS-Redfield-MultiSystem." << std::endl;
		
		// If this is the first step, write first part of header to the data file
		if(this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}
		
		// Loop through all SpinSystems to obtain SpinSpace objects
		auto systems = this->SpinSystems();
		std::vector<std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>> spaces;
		unsigned int dimensions = 0;
		for(auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			auto space = std::make_shared<SpinAPI::SpinSpace>(*(*i));
			
			// We are using hilbert space because of R, later it will be transformed into Liouville space
			space->UseSuperoperatorSpace(false);
			space->SetReactionOperatorType(this->reactionOperators);
			dimensions += space->SpaceDimensions();
			
			// Make sure to save the newly created spin space
			spaces.push_back(std::pair<std::shared_ptr<SpinAPI::SpinSystem>, std::shared_ptr<SpinAPI::SpinSpace>>( *i, space ));
		}
		
		// Now, create a matrix to hold the Liouvillian superoperator and the initial state
		arma::cx_mat L((dimensions * dimensions), (dimensions * dimensions));
		arma::cx_vec rho0((dimensions * dimensions));
		unsigned int nextDimension = 0;	// Keeps track of the dimension where the next spin space starts
		
		// Loop through the systems again to fill this matrix and vector
		for(auto i = spaces.cbegin(); i != spaces.cend(); i++)
		{
			// Make sure we have an initial state
			auto initial_states = i->first->InitialState();
			arma::cx_mat rho0HS;
			if(initial_states.size() < 1)
			{
				this->Log() << "Note: No initial state specified for spin system \"" << i->first->Name() << "\", setting the initial state to zero." << std::endl;
				rho0HS = arma::zeros<arma::cx_mat>(i->second->HilbertSpaceDimensions(), i->second->HilbertSpaceDimensions());
			}
			else
			{
				// Get the initial state for the current system
				for(auto j = initial_states.cbegin(); j != initial_states.cend(); j++)
				{
					arma::cx_mat tmp_rho0;
					if(!i->second->GetState(*j, tmp_rho0))
					{
						this->Log() << "ERROR: Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\", initial state of SpinSystem \"" << i->first->Name() << "\"." << std::endl;
						return false;
					}
					if(j == initial_states.cbegin())
						rho0HS = tmp_rho0;
					else
						rho0HS += tmp_rho0;
				}
				rho0HS /= arma::trace(rho0HS);	// The density operator should have a trace of 1
			}
			
			// Now put the initial state into the superspace vector
			arma::cx_vec rho0vec;
			if(!i->second->OperatorToSuperspace(rho0HS, rho0vec))
			{
				this->Log() << "ERROR: Failed convert initial state to superspace for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}
			
			rho0.rows(nextDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1) = rho0vec;
			
			// Next, get the Hamiltonian
			arma::cx_mat H;
			if(!i->second->Hamiltonian(H))
			{
				this->Log() << "ERROR: Failed to obtain the superspace Hamiltonian for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
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
			// CONSTRUCTING TRANSITION MATRIX "domega" OUT OF EIGENVALUES OF H0
			// ----------------------------------------------------------------
			arma::cx_mat domega(size(H));

			for (int k = 0; k < (eigen_val.size()); k++)
			{
				for (int s = 0; s < (eigen_val.size()); s++)
				{
					domega(k,s) = eigen_val(k) - eigen_val(s);
				}
			}

			
			// Redfield tensor
			arma::cx_mat R;
			R.set_size(size(kron(H,H)));
			R.zeros();
			
			if(i->first->Name() == "products")
			{

			}
			else
			{
				// ---------------------------------------------------------------
				// SETUP RELAXATION OPERATOR
				// ---------------------------------------------------------------

				// Temporary Redfield tensor
				arma::cx_mat tmp_R;
				tmp_R.set_size(size(kron(H,H)));
				tmp_R.zeros();
				
				// R Tensor array pointer for parallelization
				arma::cx_mat ** ptr_R=NULL;
				int threads;
			
				// Get number of threads for tensor pointer arrays
				#pragma omp parallel
				{
					threads = omp_get_num_threads();
				}

				ptr_R = new arma::cx_mat* [threads];
				
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

				//Correlation function set ups
				std::vector<double> tau_c_list;
				std::vector<double> ampl_list;
				arma::cx_double ampl_combined = 0.0;

				// Defining variables for parameter of user
				int terms = 0;
				int def_g = 0;
				int def_specdens = 0;
				int ops = 0;
				int coeff = 0;

				// Defining varibales for loops and storage of operators
				int num_op;
				arma::cx_mat ** ptr_Tensors=NULL;

				//Spin-Operators
				arma::cx_mat *Sz1 = new arma::cx_mat;
				arma::cx_mat *Sx1 = new arma::cx_mat;
				arma::cx_mat *Sy1 = new arma::cx_mat;
				arma::cx_mat *Sz2 = new arma::cx_mat;
				arma::cx_mat *Sx2 = new arma::cx_mat;
				arma::cx_mat *Sy2 = new arma::cx_mat;

				// T0 for rank 0 & 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - i->second->Rk1SphericalTensorXXX)
				arma::cx_mat *T0_rank_0 = new arma::cx_mat;
				arma::cx_mat *T0_rank_2 = new arma::cx_mat;

				arma::cx_mat *T0_rank_1 = new arma::cx_mat;
				arma::cx_mat *Tp_rank_1 = new arma::cx_mat;
				arma::cx_mat *Tm_rank_1 = new arma::cx_mat;

				// Tp1 & T1m for rank 2 (rank 1 neglected but can be added - see SpinSpace_operators.cpp - i->second->Rk1SphericalTensorXXX)
				arma::cx_mat *Tp1 = new arma::cx_mat;
				arma::cx_mat *Tm1 = new arma::cx_mat;

				// Tp2 & Tp2 for rank 2
				arma::cx_mat *Tp2 = new arma::cx_mat;
				arma::cx_mat *Tm2 = new arma::cx_mat;

				//arma::cx_mat ** Tensors_rotating = NULL;
				int k;
				int s;
				int l;	
		
				// ------------------------------------------------------------------
				// STARTING WITH RELAXATION MATRIX CONSTRUCTION - GOOD LUCK
				// ------------------------------------------------------------------

				this->Log() << "Starting with construction of relaxation matrix." << std::endl;
				for(auto interaction = (i->first->interactions_cbegin()); interaction != (i->first->interactions_cend()); interaction++)
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

										if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s1)->Sz() ), *s1, *Sz1))			
										{
											return false;
										}	

										if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s1)->Sx() ), *s1, *Sx1))			
										{
											return false;
										}	

										if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s1)->Sy() ), *s1, *Sy1))			
										{
											return false;
										}	
								
										// Put all tensors on pointer array
										num_op = 3;
										delete[] ptr_Tensors;
										ptr_Tensors = new arma::cx_mat * [num_op];
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
																			
										//Make field complex
										complex_field = arma::conv_to<arma::cx_vec>::from(static_field);

										// Rank 0 tensor with m=0
										if(!i->second->LRk0TensorT0((*s1), complex_field, *T0_rank_0))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=0
										if(!i->second->LRk2SphericalTensorT0((*s1), complex_field, *T0_rank_2))
										{
											this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
											continue;
										}
																					
										// Rank 2 tensor with m=1
										if(!i->second->LRk2SphericalTensorTp1((*s1), complex_field, *Tp1))
										{
											this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
											continue;
										} 

										// Rank 2 tensor with m=-1
										if(!i->second->LRk2SphericalTensorTm1((*s1), complex_field, *Tm1))
										{
											this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name()  << " and Field" << "! Skipping." << std::endl;
											continue;
										}

										// Rank 2 tensor with m=2
										if(!i->second->LRk2SphericalTensorTp2((*s1), complex_field, *Tp2))
										{
											this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and Field" << "! Skipping." << std::endl;
											continue;
										}
																			
										// Rank 2 tensor with m=-2
										if(!i->second->LRk2SphericalTensorTm2((*s1), complex_field, *Tm2))
										{
											this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and Field" << "! Skipping." << std::endl;
											continue;
										}
										
										// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
										num_op = 6;
										delete[] ptr_Tensors;
										ptr_Tensors = new arma::cx_mat* [num_op];							
										ptr_Tensors[0] = T0_rank_0;
										ptr_Tensors[1] = T0_rank_2; 
										ptr_Tensors[2] = Tm1; 
										ptr_Tensors[3] = Tp1;
										ptr_Tensors[4] = Tm2;
										ptr_Tensors[5] = Tp2;
								
										// Construct spatial spherical Tensors
										SpinAPI::Tensor inTensor(0);
										if((*interaction)->Properties()->Get("tensor",inTensor) && (*interaction)->Properties()->Get("coeff", coeff) && coeff == 1)
										{
											this->Log() << "Producing spatial shperical tensors with coupling tensor..." << std::endl;
											auto ATensor = (*interaction)->CouplingTensor();
											arma::cx_mat A(3,3); 
											arma::cx_vec Am(9);
											
											A.zeros();
											A = arma::conv_to<arma::cx_mat>::from((ATensor->LabFrame()));
											
											Am(0) =  (1.0/sqrt(3.0))	* (A(0,0) + A(1,1) + A(2,2));
											Am(1) =  (1.0/sqrt(6.0))	* (3.0 * A(2,2) - (A(0,0) + A(1,1) + A(2,2)));
											Am(2) =  0.5 				* (A(0,2) + A(2,0) + ((arma::cx_double(0.0, 1.0)) * (A(1,2) + A(2,1))));
											Am(3) = -0.5 				* (A(0,2) + A(2,0) - ((arma::cx_double(0.0, 1.0)) * (A(1,2) + A(2,1))));
											Am(4) =  0.5 				* (A(0,0) - A(1,1) - ((arma::cx_double(0.0, 1.0)) * (A(0,1) + A(1,0))));
											Am(5) =  0.5 				* (A(0,0) - A(1,1) + ((arma::cx_double(0.0, 1.0)) * (A(0,1) + A(1,0))));

											//Rank 0
											*ptr_Tensors[0] =  Am(0) 	* (*ptr_Tensors[0]);
											//Rank 2
											*ptr_Tensors[1] =  Am(1) 	* (*ptr_Tensors[1]);
											*ptr_Tensors[2] = -Am(3) 	* (*ptr_Tensors[2]); 
											*ptr_Tensors[3] = -Am(2) 	* (*ptr_Tensors[3]);
											*ptr_Tensors[4] =  Am(4) 	* (*ptr_Tensors[4]);
											*ptr_Tensors[5] =  Am(5)	* (*ptr_Tensors[5]); 
										}
										else
										{
											// Norm
											*ptr_Tensors[0] = (*ptr_Tensors[0]);
											*ptr_Tensors[1] = -(*ptr_Tensors[1]); 
											*ptr_Tensors[2] = -(*ptr_Tensors[2]); 
											*ptr_Tensors[3] = (*ptr_Tensors[3]);
											*ptr_Tensors[4] = (*ptr_Tensors[4]);
											*ptr_Tensors[5] = (*ptr_Tensors[5]);
										}
									}

									// Rotate tensors in eigenbasis of H0
									#pragma omp parallel for firstprivate(ptr_Tensors) shared(eigen_vec) num_threads(threads)
									for(k=0; k < num_op; k++)
									{
										*ptr_Tensors[k] = (eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec);
									}

									this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << std::endl;

									if((*interaction)->Properties()->Get("terms", terms ) && terms == 1)
									{
										this->Log() << "No cross-relaxation terms are requested" << std::endl;

										if((*interaction)->Properties()->Get("def_g", def_g ) && def_g == 1)
										{
											this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;

											// loop through all operators and construct R tensor in Liouville space
											#pragma omp parallel for num_threads(threads)
											for (int l=0; l < threads; l++)
											{ 
												*ptr_R[l] *= 0.0;
											}

											#pragma omp parallel for firstprivate(tmp_R, SpecDens, num_op) shared(ampl_list,tau_c_list, domega, interaction, def_specdens, ptr_R, ptr_Tensors) num_threads(threads)
											for(k=0; k < num_op;  k++)
											{
												// ----------------------------------------------------------------
												// CONSTRUCTING SPECTRAL DENSITY MATRIX
												// ----------------------------------------------------------------
												ampl_combined  = ampl_list[k] * ampl_list[k];
												
												if((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
												{
													SpecDens *= 0.0;
													
													if(!ConstructSpecDensSpecific(1,ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
														continue;
													}

													SpecDens *= (*interaction)->Prefactor();
												}
												else
												{
													SpecDens *= 0.0;
													
													if(!ConstructSpecDensSpecific(0,ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
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
												std::cout << *ptr_Tensors[k] << std::endl;
												std::cout << SpecDens << std::endl;
												if(!Redfieldtensor((*ptr_Tensors[k]),(*ptr_Tensors[k]),SpecDens,tmp_R))
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

											if((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
											{
												SpecDens *= 0.0;

												if(!ConstructSpecDensGeneral(1,ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}
											else
											{
												SpecDens *= 0.0;

												if(!ConstructSpecDensGeneral(0,ampl_list, tau_c_list, domega, SpecDens))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

												SpecDens *= (*interaction)->Prefactor();
											}
											
											#pragma omp parallel for num_threads(threads)
											for (l=0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

											#pragma omp parallel for firstprivate(tmp_R) shared(ampl_list,tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
											for(k=0; k < num_op;  k++)
											{						
												// -----------------------------------------------------------------
												// CONSTRUCTING R MATRIX 
												// -----------------------------------------------------------------

												tmp_R *= 0.0;

												if(!Redfieldtensor((*ptr_Tensors[k]),(*ptr_Tensors[k]),SpecDens,tmp_R))
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
										if((*interaction)->Properties()->Get("def_g", def_g ) && def_g == 1)
										{
											this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;
											
											#pragma omp parallel for num_threads(threads)
											for (l=0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

											// loop through all operators and construct R tensor in Liouville space	
											#pragma omp parallel for collapse(2) firstprivate(tmp_R,SpecDens, ampl_combined) shared(ampl_list,tau_c_list, domega, num_op, ptr_Tensors, ptr_R, interaction, def_specdens) num_threads(threads)
											for(k=0; k < num_op;  k++)
											{
												for(s=0; s < num_op; s++)
												{
													// ----------------------------------------------------------------
													// CONSTRUCTING SPECTRAL DENSITY MATRIX
													// ----------------------------------------------------------------

													ampl_combined  = ampl_list[k] * ampl_list[s];

													if((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
													{
														SpecDens *= 0.0;
														if(!ConstructSpecDensSpecific(1,ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}

														SpecDens *= (*interaction)->Prefactor();
													}
													else
													{
														SpecDens *= 0.0;
														
														if(!ConstructSpecDensSpecific(0,ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
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

													if(!Redfieldtensor((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;

														if(!Redfieldtensor((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
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
											// loop through all operators and construct R tensor in Liouville space	
											// ----------------------------------------------------------------
											// CONSTRUCTING SPECTRAL DENSITY MATRIX
											// ----------------------------------------------------------------
											
											this->Log() << "J is generally constructed for all operators - def_g == 0 " << std::endl;
											if((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
											{
												SpecDens *= 0.0;
											
												if(!ConstructSpecDensGeneral(1,ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}
											else
											{
												SpecDens *= 0.0;

												if(!ConstructSpecDensGeneral(0,ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}
											
											#pragma omp parallel for num_threads(threads)
											for (l=0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

											#pragma omp parallel for collapse(2) firstprivate(tmp_R) shared(ampl_list,tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
											for(k=0; k < num_op;  k++)
											{
												for(s=0; s < num_op; s++)
												{
													// ----------------------------------------------------------------
													// CONSTRUCTING R MATRIX 
													// ----------------------------------------------------------------

													tmp_R *= 0.0;
																								
													if(!Redfieldtensor((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}
													
													*ptr_R[omp_get_thread_num()] += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;

														if(!Redfieldtensor((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
														{
															this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
															continue;
														}
													}
													
													*ptr_R[omp_get_thread_num()] += tmp_R;
												}
											}	
										}
									}
							
									for (l=0; l < threads; l++)
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
											
											if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s1)->Sz() ), *s1, *Sz1))			
											{
												return false;
											}	

											if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s1)->Sx() ), *s1, *Sx1))			
											{
												return false;
											}	

											if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s1)->Sy() ), *s1, *Sy1))			
											{
												return false;
											}	
											
											if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s2)->Sz() ), *s2, *Sz2))			
											{
												return false;
											}	

											if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s2)->Sx() ), *s2, *Sx2))			
											{
												return false;
											}	

											if(!i->second->CreateOperator(arma::conv_to<arma::cx_mat>::from( (*s2)->Sy() ), *s2, *Sy2))			
											{
												return false;
											}	
											
											// Put all tensors on pointer array
											num_op = 6;
											delete[] ptr_Tensors;
											ptr_Tensors = new arma::cx_mat * [num_op];
											ptr_Tensors[0] = Sx1;
											ptr_Tensors[1] = Sy1;
											ptr_Tensors[2] = Sz1;
											ptr_Tensors[3] = Sx2;
											ptr_Tensors[4] = Sy2;
											ptr_Tensors[5] = Sz2;

											// Rotate tensors in eigenbasis of H0
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
											
											// ----------------------------------------------------------------
											// CREATION OF IRREDUCIBLE SPHERICAL TENSORS
											// ----------------------------------------------------------------			

											if(!i->second->BlRk0TensorT0(*s1, *s2, *T0_rank_0))
											{
												this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Rank 2 tensor with m=0
											if(!i->second->BlRk2SphericalTensorT0(*s1, *s2, *T0_rank_2))
											{
												this->Log() << "Failed to produce T0 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}
																						
											// Rank 2 tensor with m=1
											if(!i->second->BlRk2SphericalTensorTp1(*s1, *s2, *Tp1))
											{
												this->Log() << "Failed to produce Tp1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											} 

											// Rank 2 tensor with m=-1
											if(!i->second->BlRk2SphericalTensorTm1(*s1, *s2, *Tm1))
											{
												this->Log() << "Failed to produce Tm1 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}
														
											// Rank 2 tensor with m=2
											if(!i->second->BlRk2SphericalTensorTp2(*s1, *s2, *Tp2))
											{
												this->Log() << "Failed to produce Tp2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}
																						
											// Rank 2 tensor with m=-2
											if(!i->second->BlRk2SphericalTensorTm2(*s1, *s2, *Tm2))
											{
												this->Log() << "Failed to produce Tm2 spherical tensor between spin " << (*s1)->Name() << " and " << (*s2)->Name() << "! Skipping." << std::endl;
												continue;
											}

											// Put all tensors on pointer array and get the number of indicies for subsequent loops	- adjust number when including rank 1 tensors
											num_op = 6;				
											delete[] ptr_Tensors;									
											ptr_Tensors = new arma::cx_mat* [num_op];
											ptr_Tensors[0] = T0_rank_0;
											ptr_Tensors[1] = T0_rank_2; 
											ptr_Tensors[2] = Tm1; 
											ptr_Tensors[3] = Tp1;
											ptr_Tensors[4] = Tm2;
											ptr_Tensors[5] = Tp2;

											// Construct spatial spherical Tensors
											SpinAPI::Tensor inTensor(0);
											if((*interaction)->Properties()->Get("tensor",inTensor) && (*interaction)->Properties()->Get("coeff", coeff) && coeff == 1)
											{
												this->Log() << "Producing spatial shperical tensors with coupling tensor..." << std::endl;
												auto ATensor = (*interaction)->CouplingTensor();
												arma::cx_mat A(3,3); 
												arma::cx_vec Am(6);
												
												A.zeros();
												A = arma::conv_to<arma::cx_mat>::from((ATensor->LabFrame()));
												
												Am(0) =  (1.0/sqrt(3.0))	* (A(0,0) + A(1,1) + A(2,2));
												Am(1) =  (1.0/sqrt(6.0))	* (3.0 * A(2,2) - (A(0,0) + A(1,1) + A(2,2)));
												Am(2) =  0.5 				* (A(0,2) + A(2,0) - ((arma::cx_double(0.0, 1.0)) * (A(1,2) + A(2,1))));
												Am(3) = -0.5 				* (A(0,2) + A(2,0) + ((arma::cx_double(0.0, 1.0)) * (A(1,2) + A(2,1))));
												Am(4) =  0.5 				* (A(0,0) - A(1,1) - ((arma::cx_double(0.0, 1.0)) * (A(0,1) + A(1,0))));
												Am(5) =  0.5 				* (A(0,0) - A(1,1) + ((arma::cx_double(0.0, 1.0)) * (A(0,1) + A(1,0))));
												
												//Rank 0
												*ptr_Tensors[0] =  Am(0) 	* (*ptr_Tensors[0]);
												//Rank 2
												*ptr_Tensors[1] =  Am(1) 	* (*ptr_Tensors[1]);
												*ptr_Tensors[2] = -Am(3) 	* (*ptr_Tensors[2]); 
												*ptr_Tensors[3] = -Am(2) 	* (*ptr_Tensors[3]);
												*ptr_Tensors[4] =  Am(4) 	* (*ptr_Tensors[4]);
												*ptr_Tensors[5] =  Am(5)	* (*ptr_Tensors[5]); 
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
										for(k=0; k < num_op; k++)
										{
											*ptr_Tensors[k] = eigen_vec.t() * (*ptr_Tensors[k]) * eigen_vec;
										}	

										this->Log() << "Calculating R tensor for:" << (*interaction)->Name() << " for spin " << (*s1)->Name() << " and spin " << (*s2)->Name() << std::endl;

										if((*interaction)->Properties()->Get("def_g", def_g ) && def_g == 1)
										{
											this->Log() << "Setting J up for each operator separatly - def_g == 1 " << std::endl;

											#pragma omp parallel for num_threads(threads)
											for (l=0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}

											#pragma omp parallel for collapse(2) firstprivate(tmp_R,SpecDens,ampl_combined) shared(ampl_list,tau_c_list, domega, num_op, ptr_R, ptr_Tensors, def_specdens, interaction) num_threads(threads)
											for(k=0; k < num_op;  k++)
											{
												for(s=0; s < num_op; s++)
												{
													// ----------------------------------------------------------------
													// CONSTRUCTING SPECTRAL DENSITY MATRIX
													// ----------------------------------------------------------------

													ampl_combined  = ampl_list[k] * ampl_list[s];

													if((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
													{
														SpecDens *= 0.0;

														if(!ConstructSpecDensSpecific(1,ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}

														SpecDens *= (*interaction)->Prefactor();
													}
													else
													{
														SpecDens *= 0.0;
																
														if(!ConstructSpecDensSpecific(0,ampl_combined, static_cast<std::complex<double>>(tau_c_list[0]), domega, SpecDens))
														{
															this->Log() << "There are problems with the construction of the spectral density matrix - Please check your input." << std::endl;
															continue;
														}
																
														SpecDens *= (*interaction)->Prefactor();
													}

													ampl_combined *= 0.0;
													
													// -----------------------------------------------------------------------------------------
													// CONSTRUCTING R MATRIX 
													// -----------------------------------------------------------------------------------------
													
													tmp_R *= 0.0;

													if(!Redfieldtensor((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;
														if(!Redfieldtensor((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
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
											
											if((*interaction)->Properties()->Get("def_specdens", def_specdens) && def_specdens == 1)
											{
												SpecDens *= 0.0;

												if(!ConstructSpecDensGeneral(1,ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}
											else
											{
												SpecDens *= 0.0;

												if(!ConstructSpecDensGeneral(0,ampl_list, tau_c_list, domega, SpecDens))
												{
													this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
													continue;
												}

												SpecDens *= (*interaction)->Prefactor();
											}
											
											#pragma omp parallel for num_threads(threads)
											for (l=0; l < threads; l++)
											{
												*ptr_R[l] *= 0.0;
											}
											
											#pragma omp parallel for collapse(2) firstprivate(tmp_R) shared(ampl_list,tau_c_list, domega, SpecDens, num_op, ptr_R, ptr_Tensors) num_threads(threads)
											for(k=0; k < num_op;  k++)
											{
												for(s=0; s < num_op; s++)
												{
													// -----------------------------------------------------------------------------------------
													// CONSTRUCTING R MATRIX 
													// -----------------------------------------------------------------------------------------
													tmp_R *= 0.0;

													if(!Redfieldtensor((*ptr_Tensors[k]),(*ptr_Tensors[s]),SpecDens,tmp_R))
													{
														this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
														continue;
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;

													if (k != s)
													{
														tmp_R *= 0.0;
														if(!Redfieldtensor((*ptr_Tensors[s]),(*ptr_Tensors[k]),SpecDens,tmp_R))
														{
															this->Log() << "There are problems with the construction of the Redfield tensor - Please check your input." << std::endl;
															continue;
														}
													}

													*ptr_R[omp_get_thread_num()] += tmp_R;
												}	
											}
										}
										
										for (l=0; l < threads; l++)
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
				}
				
				for (int l = 0; l < threads; l++) 
				{
					delete ptr_R[l];
				}		
				delete [] ptr_R;
							
				for (int l = 0; l < num_op; l++) 
				{
					delete ptr_Tensors[l];
				}
				delete [] ptr_Tensors;
			}
			// ---------------------------------------------------------------
			// SETUP COMPLETE HAMILTONIAN
			// ---------------------------------------------------------------
			// Transform  H0 into superspace
			arma::cx_mat lhs;
			arma::cx_mat rhs;
			arma::cx_mat eig_val_mat;
			arma::cx_mat H_SS;

			// Constructing diagonal matrix with eigenvalues of H0
			eig_val_mat = diagmat(arma::conv_to<arma::cx_mat>::from(eigen_val));

			// Transforming into superspace
			i->second->SuperoperatorFromLeftOperator(eig_val_mat, lhs);
			i->second->SuperoperatorFromRightOperator(eig_val_mat, rhs);

		    H_SS = lhs - rhs;
			
			//Conversion into Liouville space of static Hamiltonian

			L.submat(nextDimension, nextDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1) = arma::cx_double(0.0, -1.0) * H_SS;
			
			// Then get the reaction operators
			arma::cx_mat K;
			if(!i->second->TotalReactionOperator(K))
			{
				this->Log() << "ERROR: Failed to obtain matrix representation of the reaction operators for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}

			// Rotation into eigenbasis of H0
			K = eigen_vec.t() * K * eigen_vec;

			// Transform ReactionOperator into superspace
			arma::cx_mat Klhs;
			arma::cx_mat Krhs;
			arma::cx_mat K_SS;
			
			i->second->SuperoperatorFromLeftOperator(K, Klhs);
			i->second->SuperoperatorFromRightOperator(K, Krhs);

            K_SS = Klhs + Krhs;

			//Conversion into Liouville space of reaction operator

			L.submat(nextDimension, nextDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1) -= K_SS;

			//Include R tensor in Liouville operator

			L.submat(nextDimension, nextDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1) += R;
			
			// Obtain the creation operators - note that we need to loop through the other SpinSystems again to find transitions leading into the current SpinSystem
			unsigned int nextCDimension = 0;	// Similar to nextDimension, but to keep track of first dimension for this other SpinSystem
			for(auto j = spaces.cbegin(); j != spaces.cend(); j++)
			{
				// Creation operators are off-diagonal elements
				if(j != i)
				{
					// Check all transitions whether they should produce a creation operator
					for(auto t = j->first->transitions_cbegin(); t != j->first->transitions_cend(); t++)
					{
						// Does the Transition lead into the current spin space?
						if((*t)->Target() == i->first)
						{
							// Prepare a creation operator
							arma::sp_cx_mat C;
							if(!SpinAPI::CreationOperator((*t), *(j->second), *(i->second), C, true))
							{
								this->Log() << "ERROR: Failed to obtain matrix representation of the creation operator for transition \"" << (*t)->Name() << "\"!" << std::endl;
								return false;
							}
							
							// Put it into the total Liouvillian:
							//  - The row should be that of the current spin space (the target space)
							//  - The column should be that of the source spin space (the spin system containing the Transition object)
							L.submat(nextDimension, nextCDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1, nextCDimension + (j->second->SpaceDimensions() * j->second->SpaceDimensions()) - 1) += C * (*t)->Rate();
						}
					}
				}
				
				// Move on to check next spin system for transitions into the current spin space
				nextCDimension += (j->second->SpaceDimensions() * j->second->SpaceDimensions());
			}
			
			// Move on to next spin space
			nextDimension += (i->second->SpaceDimensions() * i->second->SpaceDimensions());
		}
		
		// Write results for initial state as well (i.e. at time 0)
		this->Data() << this->RunSettings()->CurrentStep() << " 0 ";
		this->WriteStandardOutput(this->Data());
		nextDimension = 0;
		
		for(auto i = spaces.cbegin(); i != spaces.cend(); i++)
		{
			// Get the superspace result vector and convert it back to the native Hilbert space
			arma::cx_mat rho_result;
			arma::cx_vec rho_result_vec;
			rho_result_vec = rho0.rows(nextDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1);
			if(!i->second->OperatorFromSuperspace(rho_result_vec, rho_result))
			{
				this->Log() << "ERROR: Failed to convert resulting superspace-vector back to native Hilbert space for spin system \"" << i->first->Name() << "\"!" << std::endl;
				return false;
			}
		
			// Get the results
			this->GatherResults(rho_result, *(i->first), *(i->second));
	
			// Move on to next spin space
			nextDimension += (i->second->SpaceDimensions() * i->second->SpaceDimensions());
		}
		this->Data() << std::endl;
		
		// We need the propagator
		this->Log() << "Calculating the propagator..." << std::endl;
		arma::cx_mat P = arma::expmat(arma::conv_to<arma::cx_mat>::from(L) * this->timestep);
		
		// Perform the calculation
		this->Log() << "Ready to perform calculation." << std::endl;
		unsigned int steps = static_cast<unsigned int>(std::abs(this->totaltime / this->timestep));
		for(unsigned int n = 1; n <= steps; n++)
		{
			// Write first part of the data output
			this->Data() << this->RunSettings()->CurrentStep() << " ";
			this->Data() << (static_cast<double>(n) * this->timestep) << " ";
			this->WriteStandardOutput(this->Data());
			
			// Propagate (use special scope to be able to dispose of the temporary vector asap)
			{
				arma::cx_vec tmp = P * rho0;
				rho0 = tmp;
			}
			
			// Retrieve the resulting density matrix for each spin system and output the results
			nextDimension = 0;
			for(auto i = spaces.cbegin(); i != spaces.cend(); i++)
			{
				// Get the superspace result vector and convert it back to the native Hilbert space
				arma::cx_mat rho_result;
				arma::cx_vec rho_result_vec;
				rho_result_vec = rho0.rows(nextDimension, nextDimension + (i->second->SpaceDimensions() * i->second->SpaceDimensions()) - 1);
				if(!i->second->OperatorFromSuperspace(rho_result_vec, rho_result))
				{
					this->Log() << "ERROR: Failed to convert resulting superspace-vector back to native Hilbert space for spin system \"" << i->first->Name() << "\"!" << std::endl;
					return false;
				}
				
				// Get the results
				this->GatherResults(rho_result, *(i->first), *(i->second));
			
				// Move on to next spin space
				nextDimension += (i->second->SpaceDimensions() * i->second->SpaceDimensions());
			}
			
			// Terminate the line in the data file after iteration through all spin systems
			this->Data() << std::endl;
		}
		
		this->Log() << "Done with calculation." << std::endl;
		
		return true;
	}
	
	// Gathers and outputs the results from a given time-integrated density operator
	void TaskMultiStaticSSRedfieldTimeEvo::GatherResults(const arma::cx_mat& _rho, const SpinAPI::SpinSystem& _system, const SpinAPI::SpinSpace& _space)
	{
		// Loop through all states
		arma::cx_mat P;
		auto states = _system.States();
		for(auto j = states.cbegin(); j != states.cend(); j++)
		{
			if(!_space.GetState((*j), P))
			{
				this->Log() << "Failed to obtain projection matrix onto state \"" << (*j)->Name() << "\" of SpinSystem \"" << _system.Name() << "\"." << std::endl;
				continue;
			}
			
			// Return the yield for this state - note that no reaction rates are included here.
			this->Data() << std::abs(arma::trace(P * _rho)) << " ";
		}
	}
	
	// Writes the header of the data file (but can also be passed to other streams)
	void TaskMultiStaticSSRedfieldTimeEvo::WriteHeader(std::ostream& _stream)
	{
		_stream << "Step ";
		_stream << "Time(ns) ";
		this->WriteStandardOutputHeader(_stream);
		
		// Get header for each spin system
		auto systems = this->SpinSystems();
		for(auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Write each state name
			auto states = (*i)->States();
			for(auto j = states.cbegin(); j != states.cend(); j++)
				_stream << (*i)->Name() << "." << (*j)->Name() << " ";
		}
		_stream << std::endl;
	}
	
	// Validation of the required input
	bool TaskMultiStaticSSRedfieldTimeEvo::Validate()
	{
		double inputTimestep = 0.0;
		double inputTotaltime = 0.0;
		
		// Get timestep
		if(this->Properties()->Get("timestep",inputTimestep))
		{
			if(std::isfinite(inputTimestep) && inputTimestep > 0.0) {this->timestep = inputTimestep;}
			else
			{
				// We can run the calculation if an invalid timestep was specified
				return false;
			}
		}
		
		// Get totaltime
		if(this->Properties()->Get("totaltime",inputTotaltime))
		{
			if(std::isfinite(inputTotaltime) && inputTotaltime > 0.0) {this->totaltime = inputTotaltime;}
			else
			{
				// We can run the calculation if an invalid total time was specified
				return false;
			}
		}
		
		// Get the reacton operator type
		std::string str;
		if(this->Properties()->Get("reactionoperators", str))
		{
			if(str.compare("haberkorn") == 0)
			{
				this->reactionOperators = SpinAPI::ReactionOperatorType::Haberkorn;
				this->Log() << "Setting reaction operator type to Haberkorn." << std::endl;
			}
			else if(str.compare("lindblad") == 0)
			{
				this->reactionOperators = SpinAPI::ReactionOperatorType::Lindblad;
				this->Log() << "Setting reaction operator type to Lindblad." << std::endl;
			}
			else
			{
				this->Log() << "Warning: Unknown reaction operator type specified. Using default reaction operators." << std::endl;
			}
		}
		
		return true;
	}
		//Construction Refield tensor
	bool TaskMultiStaticSSRedfieldTimeEvo::Redfieldtensor(const arma::cx_mat& _op1, const arma::cx_mat& _op2, const arma::cx_mat& _specdens, arma::cx_mat& _redfieldtensor) 
	{

		arma::cx_mat one = arma::eye<arma::cx_mat>(arma::size(_specdens));

		_redfieldtensor*=0.0;

		// r = A1[a,c] * A2[d,b] * (S[c,a] + S[b,d])
		_redfieldtensor += arma::kron(((_op1)%_specdens.st()),(_op2).st());
		_redfieldtensor += arma::kron((_op1),(((_op2).st())%_specdens)); // (A2 * S.T).T = A2.T * S
	    // if b == d: r -= sum(A1[a,:] * A2[:,c] * S[c,:])
		_redfieldtensor -= arma::kron((_op1)*((_op2)%_specdens.st()),one);
		// if a == c: r -= sum(A1[d,:] * S[:,d] * A2[:,b])
		_redfieldtensor -= arma::kron(one,(((_op1)%_specdens.st())*(_op2)).st());

		return true;
	}

	bool TaskMultiStaticSSRedfieldTimeEvo::ConstructSpecDensGeneral(const int& _spectral_function, const std::vector<double>& _ampl_list, const std::vector<double>& _tau_c_list,const arma::cx_mat& _domega, arma::cx_mat& _specdens)
	{
		if(_spectral_function == 1)
		{
			// Solution  of spectral density : S = Ampl*(tau_c/(1+domega²*tauc²))
			#pragma omp for
			for (auto ii = 0; ii != _tau_c_list.size(); ii++)
			{
			   	_specdens += (static_cast<std::complex<double>>(_ampl_list[ii]) * (static_cast<std::complex<double>>(_tau_c_list[ii]) / (arma::cx_double(1.00,0.00) + (pow(_domega,2) * (pow(static_cast<std::complex<double>>(_tau_c_list[ii]),2)))))); 
			}
		}

		if(_spectral_function == 0)
		{
			// Solution of spectral density: S = Ampl/(1/tau_c - i * domega)
			#pragma omp for
			for (auto ii = 0; ii != _tau_c_list.size(); ii++)
				{
					_specdens +=  (static_cast<std::complex<double>>(_ampl_list[ii])) / ((1.00/static_cast<std::complex<double>>(_tau_c_list[ii]) - (arma::cx_double(0.0, 1.0)*_domega)));
				}	
		}

		return true;
	}

	bool TaskMultiStaticSSRedfieldTimeEvo::ConstructSpecDensSpecific(const int& _spectral_function, const std::complex<double>& _ampl, const std::complex<double>& _tau_c,const arma::cx_mat& _domega, arma::cx_mat& _specdens)
	{
		if(_spectral_function == 1)
		{
			// Solution  of spectral density : S = Ampl*(tau_c/(1+domega²*tauc²))
			_specdens = (static_cast<std::complex<double>>(_ampl) * (static_cast<std::complex<double>>(_tau_c) / (arma::cx_double(1.00,0.00) + (pow(_domega,2) * (pow(static_cast<std::complex<double>>(_tau_c),2)))))); 

		}

		if(_spectral_function == 0)
		{
			// Solution of spectral density: S = Ampl/(1/tau_c - i * domega)
			_specdens = _ampl / ((arma::cx_double(1.00, 0.00)/_tau_c) - (arma::cx_double(0.0, 1.0) * _domega));
		}

		return true;
	}
	// -----------------------------------------------------
}

