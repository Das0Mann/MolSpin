/////////////////////////////////////////////////////////////////////////
// TaskActionSpectrumHistogramRPOnlyDec implementation (RunSection module)
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// Task implemented by Siu Ying Wong and Luca Gerhards.
// (c) 2022 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "TaskActionSpectrumHistogramRPOnlyDec.h"
#include "Transition.h"
#include "Settings.h"
#include "State.h"
#include "SpinSpace.h"
#include "SpinSystem.h"
#include "Spin.h"
#include "Interaction.h"
#include "ObjectParser.h"

namespace RunSection
{
	// -----------------------------------------------------
	// TaskActionSpectrumHistogramRPOnlyDec Constructors and Destructor
	// -----------------------------------------------------
	TaskActionSpectrumHistogramRPOnlyDec::TaskActionSpectrumHistogramRPOnlyDec(const MSDParser::ObjectParser &_parser, const RunSection &_runsection) : BasicTask(_parser, _runsection)
	{
	}

	TaskActionSpectrumHistogramRPOnlyDec::~TaskActionSpectrumHistogramRPOnlyDec()
	{
	}

	// -----------------------------------------------------
	// TaskActionSpectrumHistogramRPOnlyDec protected methods
	// -----------------------------------------------------
	double conversionfactortoMHzRPOnlyDec = 5e2 / arma::datum::pi;

	bool TaskActionSpectrumHistogramRPOnlyDec::RunLocal()
	{
		this->Log() << "Running method ActionSpectrumHistogramRPOnlyDec." << std::endl;

		// If this is the first step, write first part of header to the data file
		auto bincenters = this->GetHistogramBinCenters();
		if (this->RunSettings()->CurrentStep() == 1)
		{
			this->WriteHeader(this->Data());
		}
		// ----------------------------------------------------------------
		// SETTING UP SPIN SYSTEM AND GETTING DENSITY OPERATOR
		// ----------------------------------------------------------------

		auto systems = this->SpinSystems();
		for (auto i = systems.cbegin(); i != systems.cend(); i++)
		{
			// Get a list of subspaces, make sure that we have a pair of uncoupled radicals
			auto subspaces = SpinAPI::CompleteSubspaces(*(*i));
			if (subspaces.size() < 2)
			{
				this->Log() << "Failed to obtain radicals. The spin system does not have two uncoupled subspaces! Skipping SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
				continue;
			}

			// Objects to hold the two unpaired electrons and their subspaces
			SpinAPI::spin_ptr radical[2] = {nullptr, nullptr};
			std::vector<SpinAPI::spin_ptr> subspace1;
			std::vector<SpinAPI::spin_ptr> subspace2;

			// Find two subspaces with an electron
			for (auto j = subspaces.cbegin(); j != subspaces.cend(); j++)
			{
				for (auto k = j->cbegin(); k != j->cend(); k++)
				{
					if ((*k)->Type() == SpinAPI::SpinType::Electron)
					{
						if (radical[0] == nullptr)
						{
							radical[0] = (*k);
							subspace1 = (*j);
						}
						else if (radical[1] == nullptr)
						{
							radical[1] = (*k);
							subspace2 = (*j);
						}

						break;
					}
				}
			}

			// Check whether we found a radical pair
			if (radical[0] == nullptr || radical[1] == nullptr)
			{
				this->Log() << "Failed to obtain radicals. Did not find two uncoupled electronic spins! Skipping SpinSystem \"" << (*i)->Name() << "\"." << std::endl;
				continue;
			}

			// Obtain SpinSpaces to describe the two subsystems
			SpinAPI::SpinSpace spaces[2] = {SpinAPI::SpinSpace(subspace1), SpinAPI::SpinSpace(subspace2)};
			spaces[0].UseSuperoperatorSpace(false);
			spaces[1].UseSuperoperatorSpace(false);
			spaces[0].Add((*i)->Interactions());
			spaces[1].Add((*i)->Interactions());

			// Provide information about the radicals
			this->Log() << "---------------------------------------" << std::endl;
			this->Log() << "Found radical 1 with " << subspace1.size() << " spins:" << std::endl;
			this->Log() << " - Unpaired Electron: " << radical[0]->Name() << std::endl;
			this->Log() << " - Subspace dimensions: " << spaces[0].SpaceDimensions() << std::endl;

			if (subspace1.size() > 1)
			{
				this->Log() << " - Other spins:" << std::endl;
				for (auto j = subspace1.cbegin(); j != subspace1.cend(); j++)
					if ((*j) != radical[0])
						this->Log() << "   - " << (*j)->Name() << std::endl;
			}
			else
			{
				this->Log() << " - There are no other spins." << std::endl;
			}

			this->Log() << "\nFound radical 2 with " << subspace2.size() << " spins:" << std::endl;
			this->Log() << " - Unpaired Electron: " << radical[1]->Name() << std::endl;
			this->Log() << " - Subspace dimensions: " << spaces[1].SpaceDimensions() << std::endl;

			if (subspace2.size() > 1)
			{
				this->Log() << " - Other spins:" << std::endl;
				for (auto j = subspace2.cbegin(); j != subspace2.cend(); j++)
					if ((*j) != radical[1])
						this->Log() << "   - " << (*j)->Name() << std::endl;
			}
			else
			{
				this->Log() << " - There are no other spins." << std::endl;
			}

			this->Log() << "\nSpins not included in radicals: " << ((*i)->spins_size() - subspace1.size() - subspace2.size()) << " / " << (*i)->spins_size() << std::endl;
			this->Log() << "---------------------------------------" << std::endl;
			this->Log() << "Spins that are not included in the radicals are not considered in the calculations." << std::endl;

			// Matrices to hold spin operators, eigenvectors and eigenvalue-differences
			arma::cx_mat eigenvectors[2]; // To hold eigenvectors
			arma::vec eigenvalues[2];	  // To hold eigenvalues

			std::cout << spaces[0].HilbertSpaceDimensions() << std::endl;
			std::cout << spaces[1].HilbertSpaceDimensions() << std::endl;
			arma::cx_mat Id1;
			Id1.eye(spaces[0].HilbertSpaceDimensions(), spaces[0].HilbertSpaceDimensions());
			arma::cx_mat Id2;
			Id2.eye(spaces[1].HilbertSpaceDimensions(), spaces[1].HilbertSpaceDimensions());

			arma::sp_cx_mat transitionhamiltonancomplete = arma::kron(arma::conv_to<arma::sp_cx_mat>::from(Id1), arma::conv_to<arma::sp_cx_mat>::from(Id2));
			transitionhamiltonancomplete.zeros();

			arma::cx_mat tmptransitionhamiltonian[2];

			// Fill the matrices (get Hamiltonian, diagonalize it, get spin operators...)
			for (unsigned int r = 0; r < 2; r++)
			{
				// Get the Hamiltonian of the first radical subspace
				arma::cx_mat H;

				if (!spaces[r].Hamiltonian(H))
				{
					this->Log() << "Failed to obtain Hamiltonian for radical " << r << "." << std::endl;
					continue;
				}

				// Diagonalize Hamiltonian
				if (!arma::eig_sym(eigenvalues[r], eigenvectors[r], H))
				{
					this->Log() << "Failed to diagonalize Hamiltonian for radical " << r << "." << std::endl;
					continue;
				}

				// Make transition hamiltonian with copy and pasted code.
				arma::cx_mat transitionhamiltonian = arma::zeros<arma::cx_mat>(arma::size(H));
				bool zeemanparametersalreadyfound = false;

				for (auto interaction = (*i)->interactions_cbegin(); interaction != (*i)->interactions_cend(); interaction++)
				{
					// Produce the transition Hamiltonian for each radical subsystem
					bool isZeemaninteraction = (*interaction)->Type() == SpinAPI::InteractionType::SingleSpin;
					bool containedinsubspace = spaces[r].Contains(*interaction);

					if (isZeemaninteraction && containedinsubspace)
					{
						// Inherit some properties of Zeeman interaction
						arma::vec staticfield = (*interaction)->Field();
						auto ATensor = (*interaction)->CouplingTensor();
						bool ignoretensors = (*interaction)->IgnoreTensors();
						auto spinlist = (*interaction)->Group1();

						// Get RF Field
						std::string rffielduserinput;
						arma::vec rffield(3);
						this->Properties()->Get("rf_field", rffielduserinput);

						bool useperpendicularfield = rffielduserinput == "perpendicular";

						if (useperpendicularfield)
							rffield = this->GetPerpendicular3DVectorRPOnlyDec(staticfield, true);
						else
							this->Properties()->Get("rf_field", rffield);

						// Check if a set of Zeeman interaction parameters have been found and parsed previously
						if (useperpendicularfield and zeemanparametersalreadyfound)
							throw std::logic_error("RF field vector must be supplied by user if multiple static field Zeeman interactions are provided.");
						zeemanparametersalreadyfound = true;

						// Obtain the list of spins interacting with the field, and define matrices to hold the magnetic moment operators
						arma::cx_mat Sx;
						arma::cx_mat Sy;
						arma::cx_mat Sz;

						// Loop through list of spins and append to the transition hamiltonian
						for (auto j = spinlist.cbegin(); j != spinlist.cend(); j++)
						{
							if (ignoretensors)
							{
								spaces[r].CreateOperator(arma::conv_to<arma::cx_mat>::from((*j)->Sx()), (*j), Sx);
								spaces[r].CreateOperator(arma::conv_to<arma::cx_mat>::from((*j)->Sy()), (*j), Sy);
								spaces[r].CreateOperator(arma::conv_to<arma::cx_mat>::from((*j)->Sz()), (*j), Sz);
							}
							else
							{
								spaces[r].CreateOperator(arma::conv_to<arma::cx_mat>::from((*j)->Tx()), (*j), Sx);
								spaces[r].CreateOperator(arma::conv_to<arma::cx_mat>::from((*j)->Ty()), (*j), Sy);
								spaces[r].CreateOperator(arma::conv_to<arma::cx_mat>::from((*j)->Tz()), (*j), Sz);
							}

							if (ATensor != nullptr && !IsIsotropic(*ATensor))
							{
								// Use the tensor to calculate the product S * A * B
								auto A = ATensor->LabFrame();
								transitionhamiltonian += Sx * rffield(0) * A(0, 0) + Sx * rffield(1) * A(0, 1) + Sx * rffield(2) * A(0, 2);
								transitionhamiltonian += Sy * rffield(0) * A(1, 0) + Sy * rffield(1) * A(1, 1) + Sy * rffield(2) * A(1, 2);
								transitionhamiltonian += Sz * rffield(0) * A(2, 0) + Sz * rffield(1) * A(2, 1) + Sz * rffield(2) * A(2, 2);
							}
							else
							{
								transitionhamiltonian += Sx * rffield(0) + Sy * rffield(1) + Sz * rffield(2);
							}
						}
						// Multiply by the isotropic value, if the tensor was isotropic
						if (ATensor != nullptr && IsIsotropic(*ATensor))
							transitionhamiltonian *= ATensor->Isotropic();

						// Multiply with the given prefactor (or 1.0 if none was specified)
						transitionhamiltonian *= (*interaction)->Prefactor();
						// Multiply by common prefactor (bohr magneton / hbar)
						if ((*interaction)->AddCommonPrefactor())
							transitionhamiltonian *= 8.79410005e+1;
					}
				}

				transitionhamiltonian = eigenvectors[r].t() * transitionhamiltonian * eigenvectors[r];
				// Sum up both transiton Hamiltonians
				tmptransitionhamiltonian[r] = transitionhamiltonian;
				std::cout << "tmphamiltonianloop " << tmptransitionhamiltonian[r].is_zero() << std::endl;
			}

			std::cout << "ID 1 " << Id1.is_zero() << std::endl;
			std::cout << "ID 2 " << Id2.is_zero() << std::endl;

			// ----------------------------------------------------------------
			// GET ARRAY OF ENERGY GAPS, SIMILAR TO THE "domega" matrix
			// ----------------------------------------------------------------

			transitionhamiltonancomplete = arma::kron(arma::conv_to<arma::sp_cx_mat>::from(tmptransitionhamiltonian[0]), arma::conv_to<arma::sp_cx_mat>::from(Id2)) + arma::kron(arma::conv_to<arma::sp_cx_mat>::from(Id1), arma::conv_to<arma::sp_cx_mat>::from(tmptransitionhamiltonian[1]));
			std::cout << "transitionhamiltonian " << transitionhamiltonancomplete.is_zero() << std::endl;

			// Get the all eigenvalues into one big space matrix
			arma::vec totaleigenvalues;
			totaleigenvalues = arma::kron(eigenvalues[0], arma::ones(spaces[1].HilbertSpaceDimensions())) + arma::kron(arma::ones(spaces[0].HilbertSpaceDimensions()), eigenvalues[1]);

			std::cout << "eigenvalues " << totaleigenvalues.is_zero() << std::endl;
			arma::vec energygaps = this->GetDifferencesRPOnlyDec(totaleigenvalues);
			std::cout << "energygaps " << energygaps.is_zero() << std::endl;

			// ----------------------------------------------------------------
			// EVALUATE RESONANCE EFFECTS USING EIGENVECTORS OF H0
			// ----------------------------------------------------------------

			// Produce the intial state
			arma::cx_vec rho0vec = this->GetPopulationVector(spaces[0].SpaceDimensions(), spaces[1].SpaceDimensions(), eigenvectors[0], eigenvectors[1]);
			std::cout << "rho0vec " << rho0vec.is_zero() << std::endl;

			arma::vec resonanceeffects = this->GetResonanceEffectsRPOnlyDec(rho0vec, transitionhamiltonancomplete);
			std::cout << "resonanceeffects " << resonanceeffects.is_zero() << std::endl;
			// std::cout << resonanceeffects << std::endl;
			//  ----------------------------------------------------------------
			//  BINNING PROCESS FOR HISTOGRAM
			//  ----------------------------------------------------------------
			double upperlimit;
			this->Properties()->Get("upper_limit", upperlimit);
			double binwidth;
			this->Properties()->Get("bin_width", binwidth);
			bool useMHz;
			this->Properties()->Get("units_in_MHz", useMHz);

			int numbins = arma::regspace(binwidth / 2, binwidth, upperlimit).size();
			auto normalisedheights = GetNormalisedHistogramHeightsRPOnlyDec(energygaps, resonanceeffects, binwidth, numbins, useMHz);
			std::cout << "normalisedheights " << normalisedheights.is_zero() << std::endl;
			this->Data() << this->RunSettings()->CurrentStep() << " ";

			this->Data() << normalisedheights.st(); // already contains endl
		}
		return true;
	}

	// --------------------------------------------------------------------------------------------------------------------------------------

	long long TaskActionSpectrumHistogramRPOnlyDec::GetNumDifferencesRPOnlyDec(long long dimension)
	{
		long long numgaps = dimension * (dimension - 1) / 2;
		return numgaps;
	}

	arma::vec TaskActionSpectrumHistogramRPOnlyDec::GetDifferencesRPOnlyDec(arma::vec values)
	{
		long long dimension = values.size();
		std::cout << dimension << std::endl;
		long long numgaps = GetNumDifferencesRPOnlyDec(dimension);
		std::cout << numgaps << std::endl;
		arma::vec differences(numgaps);
		int ind = 0;
		for (int k = 0; k < dimension - 1; k++)
		{
			for (int s = k + 1; s < dimension; s++)
			{
				differences(ind) = std::abs(values(k) - values(s));
				ind++;
			}
		}

		return differences;
	}

	arma::cx_vec TaskActionSpectrumHistogramRPOnlyDec::GetPopulationVector(const int &Z1, const int &Z2, const arma::cx_mat &V1, const arma::cx_mat &V2)
	{
		std::complex<double> i(0, 1);

		arma::cx_mat Sx = {{std::complex<double>(0.0, 0.0), std::complex<double>(0.5, 0.0)},
						   {std::complex<double>(0.5, 0.0), std::complex<double>(0.0, 0.0)}};

		arma::cx_mat Sy = {{std::complex<double>(0.0, 0.0), std::complex<double>(0.0, -0.5)},
						   {std::complex<double>(0.0, 0.5), std::complex<double>(0.0, 0.0)}};

		arma::cx_mat Sz = {{std::complex<double>(0.5, 0.0), std::complex<double>(0.0, 0.0)},
						   {std::complex<double>(0.0, 0.0), std::complex<double>(-0.5, 0.0)}};

		int hilbertspace = Z1 * Z2;

		// Placeholder for the resulting projection operator
		arma::cx_vec Ps(hilbertspace, arma::fill::ones);
		Ps /= 4.0;

		arma::cx_mat IdentityZ1;
		IdentityZ1.set_size(Z1, Z1);
		IdentityZ1.eye();

		arma::cx_mat IdentityZ2;
		IdentityZ2.set_size(Z2, Z2);
		IdentityZ2.eye();

		// Radical 1 - span and rotate
		arma::cx_mat kronSxZ1 = V1.t() * arma::kron(Sx, IdentityZ1) * V1;
		arma::cx_vec kronSxZ1diag = kronSxZ1.diag();
		arma::cx_mat kronSyZ1 = V1.t() * arma::kron(Sy, IdentityZ1) * V1;
		arma::cx_vec kronSyZ1diag = kronSyZ1.diag();
		arma::cx_mat kronSzZ1 = V1.t() * arma::kron(Sz, IdentityZ1) * V1;
		arma::cx_vec kronSzZ1diag = kronSzZ1.diag();

		// Radical 2 - span and rotate
		arma::cx_mat kronSxZ2 = V2.t() * arma::kron(Sx, IdentityZ2) * V2;
		arma::cx_vec kronSxZ2diag = kronSxZ2.diag();
		arma::cx_mat kronSyZ2 = V2.t() * arma::kron(Sy, IdentityZ2) * V2;
		arma::cx_vec kronSyZ2diag = kronSyZ2.diag();
		arma::cx_mat kronSzZ2 = V2.t() * arma::kron(Sz, IdentityZ2) * V2;
		arma::cx_vec kronSzZ2diag = kronSzZ2.diag();

		// Build the complete projection operator
		arma::cx_mat temp = arma::kron(arma::trans(kronSxZ1diag), kronSxZ2diag);
		Ps -= arma::conv_to<arma::cx_vec>::from(temp);
		temp = arma::kron(trans(kronSyZ1diag), kronSyZ2diag);
		Ps -= arma::conv_to<arma::cx_vec>::from(temp);
		temp = arma::kron(arma::trans(kronSzZ1diag), kronSzZ2diag);
		Ps -= arma::conv_to<arma::cx_vec>::from(temp);

		return Ps;
	}

	arma::vec TaskActionSpectrumHistogramRPOnlyDec::GetPerpendicular3DVectorRPOnlyDec(arma::vec vector, bool normalise)
	{
		arma::vec perp(3);
		perp.zeros();
		int mutuallyperpaxis = (arma::abs(vector)).index_min();
		if (mutuallyperpaxis == 0)
		{
			perp(2) = vector(1);
			perp(1) = -1 * vector(2);
		}
		if (mutuallyperpaxis == 1)
		{
			perp(2) = vector(0);
			perp(0) = -1 * vector(2);
		}
		if (mutuallyperpaxis == 2)
		{
			perp(1) = vector(0);
			perp(0) = -1 * vector(1);
		}
		if (normalise)
			perp /= arma::norm(perp);
		return perp;
	}

	arma::vec TaskActionSpectrumHistogramRPOnlyDec::GetResonanceEffectsRPOnlyDec(arma::cx_vec initialstatepopulation, arma::sp_cx_mat transitionhamiltonian)
	{
		long long dim = initialstatepopulation.n_elem;
		long long numresonances = GetNumDifferencesRPOnlyDec(dim);
		arma::vec resonanceeffects(numresonances);

		// Fill the array of resonance effects one element at the time
		int ind_resonance = 0;
		for (int k = 0; k < dim - 1; k++)
		{
			for (int s = k + 1; s < dim; s++)
			{
				double sqrtoftransprob = std::abs(static_cast<std::complex<double>>(transitionhamiltonian(k, s)));
				double value = std::abs(initialstatepopulation(k) - initialstatepopulation(s)) * (sqrtoftransprob * sqrtoftransprob);
				resonanceeffects(ind_resonance) = value;
				ind_resonance++;
			}
		}
		return resonanceeffects;
	}

	double TaskActionSpectrumHistogramRPOnlyDec::GetHistogramPrefactorRPOnlyDec(bool useMHz_units)
	{
		double prefactor;
		if (useMHz_units)
			prefactor = conversionfactortoMHzRPOnlyDec;
		else
			prefactor = 1e9;
		return prefactor;
	}

	arma::vec TaskActionSpectrumHistogramRPOnlyDec::GetNormalisedHistogramHeightsRPOnlyDec(arma::vec energygaps, arma::vec resonanceeffects, double binwidth, int numbins, bool useMHz_units)
	{
		double prefactor = GetHistogramPrefactorRPOnlyDec(useMHz_units);
		long long numresonances = energygaps.size();
		arma::vec heights(numbins);
		heights.zeros();
		for (int s = 0; s < numresonances; s++)
		{
			double gap = std::abs(energygaps(s)) * prefactor;
			double effect = std::abs(resonanceeffects(s));
			auto binnumber = floor(gap / binwidth); // only works for uniform binwidth
			heights(binnumber) += effect;
		}
		return heights / arma::sum(heights);
	}

	// Writes the header of the data file (but can also be passed to other streams)
	void TaskActionSpectrumHistogramRPOnlyDec::TaskActionSpectrumHistogramRPOnlyDec::WriteHeader(std::ostream &_stream)
	{
		_stream << "Step ";
		this->WriteStandardOutputHeader(_stream);
		auto bincenters = TaskActionSpectrumHistogramRPOnlyDec::GetHistogramBinCenters();
		_stream << bincenters.st(); // already contains endl
	}

	// Validation
	bool TaskActionSpectrumHistogramRPOnlyDec::Validate()
	{

		return true;
	}

	arma::vec TaskActionSpectrumHistogramRPOnlyDec::GetHistogramBinCenters()
	{
		double upperlimit;
		this->Properties()->Get("upper_limit", upperlimit);
		double binwidth;
		this->Properties()->Get("bin_width", binwidth);
		auto bincenters = arma::regspace(binwidth / 2, binwidth, upperlimit);
		return bincenters;
	}
}
