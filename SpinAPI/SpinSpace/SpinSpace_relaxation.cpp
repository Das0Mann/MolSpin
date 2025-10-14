/////////////////////////////////////////////////////////////////////////
// SpinSpace class (SpinAPI Module)
// ------------------
// This source file contains methods related to relaxation operators,
// which are described by Operator objects.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
namespace SpinAPI
{
	// -----------------------------------------------------
	// Relaxation operators
	// -----------------------------------------------------
	bool SpinSpace::RelaxationOperator(const operator_ptr &_operator, arma::cx_mat &_out) const
	{
		// Make sure that we have a valid operator object
		if (_operator == nullptr || !_operator->IsValid())
			return false;

		// The Lindblad relaxation operators can only be constructed in the superspace
		if (!this->useSuperspace)
			return false;

		if (_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();

			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;																			  // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;																			  // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate2();

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate3();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationLindbladDoubleSpin)
		{
			auto spins = _operator->Spins();

			// Vectors to hold the spin operators for each spin
			std::vector<arma::cx_mat> Sp_operators;
			std::vector<arma::cx_mat> Sm_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::cx_mat Sptmp, Smtmp;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sp()), (*i), Sptmp);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sm()), (*i), Smtmp);

				// Store the operators in the vectors
				Sp_operators.push_back(Sptmp);
				Sm_operators.push_back(Smtmp);
			}

			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;																			  // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;																			  // The right-hand operator (. * conj(P)*P)

			for (auto it1 = spins.cbegin(); it1 != spins.cend(); ++it1)
			{
				for (auto it2 = spins.cbegin(); it2 != spins.cend(); ++it2)
				{
					if (it1 == it2)
						continue; // Skip self-relaxation terms

					// Get the actual spin objects
					auto spin1 = *it1;
					auto spin2 = *it2;

					// Define cross-relaxation operators
					arma::cx_mat L_minus = Sm_operators[it1 - spins.cbegin()] * Sp_operators[it2 - spins.cbegin()]; // S- I+
					arma::cx_mat L_plus = Sp_operators[it1 - spins.cbegin()] * Sm_operators[it2 - spins.cbegin()];	// S+ I-

					// Apply Lindblad terms for L_minus (cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_minus, L_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_minus.t() * L_minus, PL) || !this->SuperoperatorFromRightOperator(L_minus.t() * L_minus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

					// Apply Lindblad terms for L_plus (reverse cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_plus, L_plus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_plus.t() * L_plus, PL) || !this->SuperoperatorFromRightOperator(L_plus.t() * L_plus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
				}
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationDephasing)
		{
			auto spins = _operator->Spins();

			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator

			arma::cx_mat PRsinglet; // The right-hand operator (. * conj(P)*P)
			arma::cx_mat PLsinglet; // The left-hand operator (conj(P)*P * .)

			arma::cx_mat PRtriplet; // The right-hand operator (. * conj(P)*P)
			arma::cx_mat PLtriplet; // The left-hand operator (conj(P)*P * .)

			arma::cx_mat Psinglet;
			arma::cx_mat Ptriplet;

			// Vectors to hold the spin operators for each spin
			std::vector<arma::cx_mat> Sx_operators;
			std::vector<arma::cx_mat> Sy_operators;
			std::vector<arma::cx_mat> Sz_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::cx_mat Sxtmp, Sytmp, Sztmp;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sxtmp);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sytmp);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sztmp);

				// Store the operators in the vectors
				Sx_operators.push_back(Sxtmp);
				Sy_operators.push_back(Sytmp);
				Sz_operators.push_back(Sztmp);
			}

			arma::cx_mat E;
			E.set_size(size(Sx_operators[0]));
			E.eye();

			Psinglet = (1.0 / 4.0) * E - (Sx_operators[0] * Sx_operators[1] + Sy_operators[0] * Sy_operators[1] + Sz_operators[0] * Sz_operators[1]);
			Ptriplet = E - Psinglet;

			if (!this->SuperoperatorFromLeftOperator(Psinglet, PLsinglet) || !this->SuperoperatorFromLeftOperator(Ptriplet, PLtriplet) || !this->SuperoperatorFromRightOperator(Psinglet, PRsinglet) || !this->SuperoperatorFromRightOperator(Ptriplet, PRtriplet))
				return false;

			P = -1.0 * _operator->Rate1() * (PLsinglet * PRtriplet + PLtriplet * PRsinglet);

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationRandomFields)
		{
			// Adapted from Kattnig et al. (2016) DOI: 10.1088/1367-2630/18/6/063007

			auto spins = _operator->Spins();

			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat E = arma::eye<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());	  // Unity operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB))
					return false;
				P += -1.0 * _operator->Rate1() * ((3 / 2) * E - PB);

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB))
					return false;
				P += -1.0 * _operator->Rate2() * ((3 / 2) * E - PB);

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB))
					return false;
				P += -1.0 * _operator->Rate3() * ((3 / 2) * E - PB);
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT1)
		{
			// Routine to calculate T1 relaxation (longitudinal)
			auto spins = _operator->Spins();

			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat S_plus;
			arma::cx_mat S_minus;
			arma::cx_mat PB; // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL; // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR; // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sp()), (*i), S_plus);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sm()), (*i), S_minus);

				// Get the contribution from T1 relaxation
				if (!this->SuperoperatorFromOperators(S_plus, S_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(S_minus.t() * S_plus, PL) || !this->SuperoperatorFromRightOperator(S_minus.t() * S_plus, PR))
					return false;

				// \mathcal{L}_{T1}(\rho) = \frac{1}{T1} (2S_+ \rho S_- - \{S_-S_+, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT2)
		{
			// Routine to calculate T2 relaxation (transverse)
			auto spins = _operator->Spins();

			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;																			  // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;																			  // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

				// Get the contribution from T2 relaxation
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				// \mathcal{L}_{T2}(\rho) = \frac{1}{T2} (2S_z \rho S_z - \{S_z^2, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::Unspecified)
		{
			// Cannot construct operator if the type is not specified
			return false;
		}
		else
		{
			std::cout << "Cannot construct relaxation operator for Operator " << _operator->Name() << "!" << std::endl;
			return false;
		}

		return true;
	}

	// Sparse version of the RelaxationOperator method
	bool SpinSpace::RelaxationOperator(const operator_ptr &_operator, arma::sp_cx_mat &_out) const
	{
		// Make sure that we have a valid operator object
		if (_operator == nullptr || !_operator->IsValid())
			return false;

		// The Lindblad relaxation operators can only be constructed in the superspace
		if (!this->useSuperspace)
			return false;

		if (_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																	   // Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																	   // The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																	   // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate2();

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate3();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationLindbladDoubleSpin)
		{
			auto spins = _operator->Spins();

			// Vectors to hold the spin operators for each spin
			std::vector<arma::sp_cx_mat> Sp_operators;
			std::vector<arma::sp_cx_mat> Sm_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::sp_cx_mat Sptmp, Smtmp;

				// Create the spin operators
				this->CreateOperator((*i)->Sp(), (*i), Sptmp);
				this->CreateOperator((*i)->Sm(), (*i), Smtmp);

				// Store the operators in the vectors
				Sp_operators.push_back(Sptmp);
				Sm_operators.push_back(Smtmp);
			}

			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																					// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																					// The right-hand operator (. * conj(P)*P)

			for (auto it1 = spins.cbegin(); it1 != spins.cend(); ++it1)
			{
				for (auto it2 = spins.cbegin(); it2 != spins.cend(); ++it2)
				{
					if (it1 == it2)
						continue; // Skip self-relaxation terms

					// Get the actual spin objects
					auto spin1 = *it1;
					auto spin2 = *it2;

					// Define cross-relaxation operators
					arma::sp_cx_mat L_minus = Sm_operators[it1 - spins.cbegin()] * Sp_operators[it2 - spins.cbegin()]; // S- I+
					arma::sp_cx_mat L_plus = Sp_operators[it1 - spins.cbegin()] * Sm_operators[it2 - spins.cbegin()];  // S+ I-

					// Apply Lindblad terms for L_minus (cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_minus, L_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_minus.t() * L_minus, PL) || !this->SuperoperatorFromRightOperator(L_minus.t() * L_minus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

					// Apply Lindblad terms for L_plus (reverse cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_plus, L_plus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_plus.t() * L_plus, PL) || !this->SuperoperatorFromRightOperator(L_plus.t() * L_plus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
				}
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationDephasing)
		{

			auto spins = _operator->Spins();

			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator

			arma::sp_cx_mat PRsinglet; // The right-hand operator (. * conj(P)*P)
			arma::sp_cx_mat PLsinglet; // The left-hand operator (conj(P)*P * .)

			arma::sp_cx_mat PRtriplet; // The right-hand operator (. * conj(P)*P)
			arma::sp_cx_mat PLtriplet; // The left-hand operator (conj(P)*P * .)

			arma::sp_cx_mat Psinglet;
			arma::sp_cx_mat Ptriplet;

			// Vectors to hold the spin operators for each spin
			std::vector<arma::sp_cx_mat> Sx_operators;
			std::vector<arma::sp_cx_mat> Sy_operators;
			std::vector<arma::sp_cx_mat> Sz_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::sp_cx_mat Sxtmp, Sytmp, Sztmp;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sxtmp);
				this->CreateOperator((*i)->Sy(), (*i), Sytmp);
				this->CreateOperator((*i)->Sz(), (*i), Sztmp);

				// Store the operators in the vectors
				Sx_operators.push_back(Sxtmp);
				Sy_operators.push_back(Sytmp);
				Sz_operators.push_back(Sztmp);
			}

			arma::sp_cx_mat E;
			E.set_size(size(Sx_operators[0]));
			E.eye();

			Psinglet = (1.0 / 4.0) * E - (Sx_operators[0] * Sx_operators[1] + Sy_operators[0] * Sy_operators[1] + Sz_operators[0] * Sz_operators[1]);
			Ptriplet = E - Psinglet;

			if (!this->SuperoperatorFromLeftOperator(Psinglet, PLsinglet) || !this->SuperoperatorFromLeftOperator(Ptriplet, PLtriplet) || !this->SuperoperatorFromRightOperator(Psinglet, PRsinglet) || !this->SuperoperatorFromRightOperator(Ptriplet, PRtriplet))
				return false;

			P = -1.0 * _operator->Rate1() * (PLsinglet * PRtriplet + PLtriplet * PRsinglet);

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationRandomFields)
		{
			// Adapted from Kattnig et al. (2016) DOI: 10.1088/1367-2630/18/6/063007

			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat E = arma::eye<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());	// Unity operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB))
					return false;
				P += -1.0 * _operator->Rate1() * ((3 / 2) * E - PB);

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB))
					return false;
				P += -1.0 * _operator->Rate2() * ((3 / 2) * E - PB);

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB))
					return false;
				P += -1.0 * _operator->Rate3() * ((3 / 2) * E - PB);
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT1)
		{
			// Routine to calculate T1 relaxation (longitudinal)
			auto spins = _operator->Spins();

			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat S_plus;
			arma::sp_cx_mat S_minus;
			arma::sp_cx_mat PB; // Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL; // The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR; // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sp(), (*i), S_plus);
				this->CreateOperator((*i)->Sm(), (*i), S_minus);

				// Get the contribution from T1 relaxation
				if (!this->SuperoperatorFromOperators(S_plus, S_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(S_minus.t() * S_plus, PL) || !this->SuperoperatorFromRightOperator(S_minus.t() * S_plus, PR))
					return false;

				// \mathcal{L}_{T1}(\rho) = \frac{1}{T1} (2S_+ \rho S_- - \{S_-S_+, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT2)
		{
			// Routine to calculate T2 relaxation (transverse)
			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																					// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																					// The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Get the contribution from T2 relaxation
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				// \mathcal{L}_{T2}(\rho) = \frac{1}{T2} (2S_z \rho S_z - \{S_z^2, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::Unspecified)
		{
			// Cannot construct operator if the type is not specified
			return false;
		}
		else
		{
			std::cout << "Cannot construct relaxation operator for Operator " << _operator->Name() << "!" << std::endl;
			return false;
		}

		return true;
	}

	// --------------------------------------------------------------
	// Relaxation operators when a unitary transformation is required
	// --------------------------------------------------------------

	bool SpinSpace::RelaxationOperatorFrameChange(const operator_ptr &_operator, arma::cx_mat _rotationmatrix, arma::cx_mat &_out) const
	{
		// Make sure that we have a valid operator object
		if (_operator == nullptr || !_operator->IsValid())
			return false;

		// The Lindblad relaxation operators can only be constructed in the superspace
		if (!this->useSuperspace)
			return false;

		if (_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();

			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;																			  // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;																			  // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

				// Rotate into the new frame with given unitary rotation matrix
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate2();

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate3();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationLindbladDoubleSpin)
		{
			auto spins = _operator->Spins();

			// Vectors to hold the spin operators for each spin
			std::vector<arma::cx_mat> Sp_operators;
			std::vector<arma::cx_mat> Sm_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::cx_mat Sptmp, Smtmp;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sp()), (*i), Sptmp);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sm()), (*i), Smtmp);

				Sptmp = _rotationmatrix.t() * Sptmp * _rotationmatrix;
				Smtmp = _rotationmatrix.t() * Smtmp * _rotationmatrix;

				// Store the operators in the vectors
				Sp_operators.push_back(Sptmp);
				Sm_operators.push_back(Smtmp);
			}

			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;																			  // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;																			  // The right-hand operator (. * conj(P)*P)

			for (auto it1 = spins.cbegin(); it1 != spins.cend(); ++it1)
			{
				for (auto it2 = spins.cbegin(); it2 != spins.cend(); ++it2)
				{
					if (it1 == it2)
						continue; // Skip self-relaxation terms

					// Get the actual spin objects
					auto spin1 = *it1;
					auto spin2 = *it2;

					// Define cross-relaxation operators
					arma::cx_mat L_minus = Sm_operators[it1 - spins.cbegin()] * Sp_operators[it2 - spins.cbegin()]; // S- I+
					arma::cx_mat L_plus = Sp_operators[it1 - spins.cbegin()] * Sm_operators[it2 - spins.cbegin()];	// S+ I-

					// Apply Lindblad terms for L_minus (cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_minus, L_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_minus.t() * L_minus, PL) || !this->SuperoperatorFromRightOperator(L_minus.t() * L_minus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

					// Apply Lindblad terms for L_plus (reverse cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_plus, L_plus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_plus.t() * L_plus, PL) || !this->SuperoperatorFromRightOperator(L_plus.t() * L_plus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
				}
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationDephasing)
		{
			auto spins = _operator->Spins();

			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator

			arma::cx_mat PRsinglet; // The right-hand operator (. * conj(P)*P)
			arma::cx_mat PLsinglet; // The left-hand operator (conj(P)*P * .)

			arma::cx_mat PRtriplet; // The right-hand operator (. * conj(P)*P)
			arma::cx_mat PLtriplet; // The left-hand operator (conj(P)*P * .)

			arma::cx_mat Psinglet;
			arma::cx_mat Ptriplet;

			// Vectors to hold the spin operators for each spin
			std::vector<arma::cx_mat> Sx_operators;
			std::vector<arma::cx_mat> Sy_operators;
			std::vector<arma::cx_mat> Sz_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::cx_mat Sxtmp, Sytmp, Sztmp;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sxtmp);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sytmp);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sztmp);

				// Store the operators in the vectors
				Sx_operators.push_back(Sxtmp);
				Sy_operators.push_back(Sytmp);
				Sz_operators.push_back(Sztmp);
			}

			arma::cx_mat E;
			E.set_size(arma::size(Sx_operators[0]));
			E.eye();

			Psinglet = (1.0 / 4.0) * E - (Sx_operators[0] * Sx_operators[1] + Sy_operators[0] * Sy_operators[1] + Sz_operators[0] * Sz_operators[1]);
			Ptriplet = E - Psinglet;

			// Rotate into the new frame with given unitary rotation matrix
			Psinglet = _rotationmatrix.t() * Psinglet * _rotationmatrix;
			Ptriplet = _rotationmatrix.t() * Ptriplet * _rotationmatrix;

			if (!this->SuperoperatorFromLeftOperator(Psinglet, PLsinglet) || !this->SuperoperatorFromLeftOperator(Ptriplet, PLtriplet) || !this->SuperoperatorFromRightOperator(Psinglet, PRsinglet) || !this->SuperoperatorFromRightOperator(Ptriplet, PRtriplet))
				return false;

			P = -1.0 * _operator->Rate1() * (PLsinglet * PRtriplet + PLtriplet * PRsinglet);

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationRandomFields)
		{
			// Adapted from Kattnig et al. (2016) DOI: 10.1088/1367-2630/18/6/063007

			auto spins = _operator->Spins();

			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat E = arma::eye<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());	  // Unity operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

				// Rotate into the eigenbasis of the spin system
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB))
					return false;
				P += -1.0 * _operator->Rate1() * ((3 / 2) * E - PB);

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB))
					return false;
				P += -1.0 * _operator->Rate2() * ((3 / 2) * E - PB);

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB))
					return false;
				P += -1.0 * _operator->Rate3() * ((3 / 2) * E - PB);
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT1)
		{
			// Routine to calculate T1 relaxation (longitudinal)
			auto spins = _operator->Spins();

			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat S_plus;
			arma::cx_mat S_minus;
			arma::cx_mat PB; // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL; // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR; // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sp()), (*i), S_plus);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sm()), (*i), S_minus);

				// Rotate into the eigenbasis of the spin system
				S_plus = _rotationmatrix.t() * S_plus * _rotationmatrix;
				S_minus = _rotationmatrix.t() * S_minus * _rotationmatrix;

				// Get the contribution from T1 relaxation
				if (!this->SuperoperatorFromOperators(S_plus, S_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(S_minus.t() * S_plus, PL) || !this->SuperoperatorFromRightOperator(S_minus.t() * S_plus, PR))
					return false;

				// \mathcal{L}_{T1}(\rho) = \frac{1}{T1} (2S_+ \rho S_- - \{S_-S_+, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT2)
		{
			// Routine to calculate T2 relaxation (transverse)
			auto spins = _operator->Spins();

			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::cx_mat PB;																			  // Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;																			  // The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;																			  // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);

				// Rotate into the eigenbasis of the spin system
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from T2 relaxation
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				// \mathcal{L}_{T2}(\rho) = \frac{1}{T2} (2S_z \rho S_z - \{S_z^2, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::Unspecified)
		{
			// Cannot construct operator if the type is not specified
			return false;
		}
		else
		{
			std::cout << "Cannot construct relaxation operator for Operator " << _operator->Name() << "!" << std::endl;
			return false;
		}

		return true;
	}

	// Sparse version of the RelaxationOperatorFrameChange method
	bool SpinSpace::RelaxationOperatorFrameChange(const operator_ptr &_operator, arma::cx_mat _rotationmatrix, arma::sp_cx_mat &_out) const
	{
		// Make sure that we have a valid operator object
		if (_operator == nullptr || !_operator->IsValid())
			return false;

		// The Lindblad relaxation operators can only be constructed in the superspace
		if (!this->useSuperspace)
			return false;

		if (_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																	   // Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																	   // The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																	   // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Rotate into the new frame with given unitary rotation matrix
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate2();

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate3();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationLindbladDoubleSpin)
		{
			auto spins = _operator->Spins();

			// Vectors to hold the spin operators for each spin
			std::vector<arma::sp_cx_mat> Sp_operators;
			std::vector<arma::sp_cx_mat> Sm_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::sp_cx_mat Sptmp, Smtmp;

				// Create the spin operators
				this->CreateOperator((*i)->Sp(), (*i), Sptmp);
				this->CreateOperator((*i)->Sm(), (*i), Smtmp);

				Sptmp = _rotationmatrix.t() * Sptmp * _rotationmatrix;
				Smtmp = _rotationmatrix.t() * Smtmp * _rotationmatrix;

				// Store the operators in the vectors
				Sp_operators.push_back(Sptmp);
				Sm_operators.push_back(Smtmp);
			}

			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																					// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																					// The right-hand operator (. * conj(P)*P)

			for (auto it1 = spins.cbegin(); it1 != spins.cend(); ++it1)
			{
				for (auto it2 = spins.cbegin(); it2 != spins.cend(); ++it2)
				{
					if (it1 == it2)
						continue; // Skip self-relaxation terms

					// Get the actual spin objects
					auto spin1 = *it1;
					auto spin2 = *it2;

					// Define cross-relaxation operators
					arma::sp_cx_mat L_minus = Sm_operators[it1 - spins.cbegin()] * Sp_operators[it2 - spins.cbegin()]; // S- I+
					arma::sp_cx_mat L_plus = Sp_operators[it1 - spins.cbegin()] * Sm_operators[it2 - spins.cbegin()];  // S+ I-

					// Apply Lindblad terms for L_minus (cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_minus, L_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_minus.t() * L_minus, PL) || !this->SuperoperatorFromRightOperator(L_minus.t() * L_minus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

					// Apply Lindblad terms for L_plus (reverse cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_plus, L_plus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_plus.t() * L_plus, PL) || !this->SuperoperatorFromRightOperator(L_plus.t() * L_plus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
				}
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationDephasing)
		{

			auto spins = _operator->Spins();

			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator

			arma::sp_cx_mat PRsinglet; // The right-hand operator (. * conj(P)*P)
			arma::sp_cx_mat PLsinglet; // The left-hand operator (conj(P)*P * .)

			arma::sp_cx_mat PRtriplet; // The right-hand operator (. * conj(P)*P)
			arma::sp_cx_mat PLtriplet; // The left-hand operator (conj(P)*P * .)

			arma::sp_cx_mat Psinglet;
			arma::sp_cx_mat Ptriplet;

			// Vectors to hold the spin operators for each spin
			std::vector<arma::sp_cx_mat> Sx_operators;
			std::vector<arma::sp_cx_mat> Sy_operators;
			std::vector<arma::sp_cx_mat> Sz_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::sp_cx_mat Sxtmp, Sytmp, Sztmp;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sxtmp);
				this->CreateOperator((*i)->Sy(), (*i), Sytmp);
				this->CreateOperator((*i)->Sz(), (*i), Sztmp);

				// Store the operators in the vectors
				Sx_operators.push_back(Sxtmp);
				Sy_operators.push_back(Sytmp);
				Sz_operators.push_back(Sztmp);
			}

			arma::sp_cx_mat E;
			E.set_size(arma::size(Sx_operators[0]));
			E.eye();

			Psinglet = (1.0 / 4.0) * E - (Sx_operators[0] * Sx_operators[1] + Sy_operators[0] * Sy_operators[1] + Sz_operators[0] * Sz_operators[1]);
			Ptriplet = E - Psinglet;

			// Rotate into the new frame with given unitary rotation matrix
			Psinglet = _rotationmatrix.t() * Psinglet * _rotationmatrix;
			Ptriplet = _rotationmatrix.t() * Ptriplet * _rotationmatrix;

			if (!this->SuperoperatorFromLeftOperator(Psinglet, PLsinglet) || !this->SuperoperatorFromLeftOperator(Ptriplet, PLtriplet) || !this->SuperoperatorFromRightOperator(Psinglet, PRsinglet) || !this->SuperoperatorFromRightOperator(Ptriplet, PRtriplet))
				return false;

			P = -1.0 * _operator->Rate1() * (PLsinglet * PRtriplet + PLtriplet * PRsinglet);

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationRandomFields)
		{
			// Adapted from Kattnig et al. (2016) DOI: 10.1088/1367-2630/18/6/063007

			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat E = arma::eye<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());	// Unity operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Rotate into the eigenbasis of the spin system
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB))
					return false;
				P += -1.0 * _operator->Rate1() * ((3 / 2) * E - PB);

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB))
					return false;
				P += -1.0 * _operator->Rate2() * ((3 / 2) * E - PB);

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB))
					return false;
				P += -1.0 * _operator->Rate3() * ((3 / 2) * E - PB);
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT1)
		{
			// Routine to calculate T1 relaxation (longitudinal)
			auto spins = _operator->Spins();

			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat S_plus;
			arma::sp_cx_mat S_minus;
			arma::sp_cx_mat PB; // Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL; // The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR; // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sp(), (*i), S_plus);
				this->CreateOperator((*i)->Sm(), (*i), S_minus);

				// Rotate into the eigenbasis of the spin system
				S_plus = _rotationmatrix.t() * S_plus * _rotationmatrix;
				S_minus = _rotationmatrix.t() * S_minus * _rotationmatrix;

				// Get the contribution from T1 relaxation
				if (!this->SuperoperatorFromOperators(S_plus, S_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(S_minus.t() * S_plus, PL) || !this->SuperoperatorFromRightOperator(S_minus.t() * S_plus, PR))
					return false;

				// \mathcal{L}_{T1}(\rho) = \frac{1}{T1} (2S_+ \rho S_- - \{S_-S_+, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT2)
		{
			// Routine to calculate T2 relaxation (transverse)
			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																					// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																					// The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Rotate into the eigenbasis of the spin system
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from T2 relaxation
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				// \mathcal{L}_{T2}(\rho) = \frac{1}{T2} (2S_z \rho S_z - \{S_z^2, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::Unspecified)
		{
			// Cannot construct operator if the type is not specified
			return false;
		}
		else
		{
			std::cout << "Cannot construct relaxation operator for Operator " << _operator->Name() << "!" << std::endl;
			return false;
		}

		return true;
	}

	// Fully sparse version of the RelaxationOperatorFrameChange method
	bool SpinSpace::RelaxationOperatorFrameChange(const operator_ptr &_operator, arma::sp_cx_mat _rotationmatrix, arma::sp_cx_mat &_out) const
	{
		// Make sure that we have a valid operator object
		if (_operator == nullptr || !_operator->IsValid())
			return false;

		// The Lindblad relaxation operators can only be constructed in the superspace
		if (!this->useSuperspace)
			return false;

		if (_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																	   // Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																	   // The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																	   // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Rotate into the new frame with given unitary rotation matrix
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate2();

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR) / 2.0) * _operator->Rate3();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationLindbladDoubleSpin)
		{
			auto spins = _operator->Spins();

			// Vectors to hold the spin operators for each spin
			std::vector<arma::sp_cx_mat> Sp_operators;
			std::vector<arma::sp_cx_mat> Sm_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::sp_cx_mat Sptmp, Smtmp;

				// Create the spin operators
				this->CreateOperator((*i)->Sp(), (*i), Sptmp);
				this->CreateOperator((*i)->Sm(), (*i), Smtmp);

				Sptmp = _rotationmatrix.t() * Sptmp * _rotationmatrix;
				Smtmp = _rotationmatrix.t() * Smtmp * _rotationmatrix;

				// Store the operators in the vectors
				Sp_operators.push_back(Sptmp);
				Sm_operators.push_back(Smtmp);
			}

			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																					// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																					// The right-hand operator (. * conj(P)*P)

			for (auto it1 = spins.cbegin(); it1 != spins.cend(); ++it1)
			{
				for (auto it2 = spins.cbegin(); it2 != spins.cend(); ++it2)
				{
					if (it1 == it2)
						continue; // Skip self-relaxation terms

					// Get the actual spin objects
					auto spin1 = *it1;
					auto spin2 = *it2;

					// Define cross-relaxation operators
					arma::sp_cx_mat L_minus = Sm_operators[it1 - spins.cbegin()] * Sp_operators[it2 - spins.cbegin()]; // S- I+
					arma::sp_cx_mat L_plus = Sp_operators[it1 - spins.cbegin()] * Sm_operators[it2 - spins.cbegin()];  // S+ I-

					// Apply Lindblad terms for L_minus (cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_minus, L_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_minus.t() * L_minus, PL) || !this->SuperoperatorFromRightOperator(L_minus.t() * L_minus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();

					// Apply Lindblad terms for L_plus (reverse cross-relaxation)
					if (!this->SuperoperatorFromOperators(L_plus, L_plus.t(), PB) || !this->SuperoperatorFromLeftOperator(L_plus.t() * L_plus, PL) || !this->SuperoperatorFromRightOperator(L_plus.t() * L_plus, PR))
						return false;

					P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
				}
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationDephasing)
		{

			auto spins = _operator->Spins();

			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator

			arma::sp_cx_mat PRsinglet; // The right-hand operator (. * conj(P)*P)
			arma::sp_cx_mat PLsinglet; // The left-hand operator (conj(P)*P * .)

			arma::sp_cx_mat PRtriplet; // The right-hand operator (. * conj(P)*P)
			arma::sp_cx_mat PLtriplet; // The left-hand operator (conj(P)*P * .)

			arma::sp_cx_mat Psinglet;
			arma::sp_cx_mat Ptriplet;

			// Vectors to hold the spin operators for each spin
			std::vector<arma::sp_cx_mat> Sx_operators;
			std::vector<arma::sp_cx_mat> Sy_operators;
			std::vector<arma::sp_cx_mat> Sz_operators;

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Temporary matrices to hold the operators
				arma::sp_cx_mat Sxtmp, Sytmp, Sztmp;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sxtmp);
				this->CreateOperator((*i)->Sy(), (*i), Sytmp);
				this->CreateOperator((*i)->Sz(), (*i), Sztmp);

				// Store the operators in the vectors
				Sx_operators.push_back(Sxtmp);
				Sy_operators.push_back(Sytmp);
				Sz_operators.push_back(Sztmp);
			}

			arma::sp_cx_mat E;
			E.set_size(arma::size(Sx_operators[0]));
			E.eye();

			Psinglet = (1.0 / 4.0) * E - (Sx_operators[0] * Sx_operators[1] + Sy_operators[0] * Sy_operators[1] + Sz_operators[0] * Sz_operators[1]);
			Ptriplet = E - Psinglet;

			// Rotate into the new frame with given unitary rotation matrix
			Psinglet = _rotationmatrix.t() * Psinglet * _rotationmatrix;
			Ptriplet = _rotationmatrix.t() * Ptriplet * _rotationmatrix;

			if (!this->SuperoperatorFromLeftOperator(Psinglet, PLsinglet) || !this->SuperoperatorFromLeftOperator(Ptriplet, PLtriplet) || !this->SuperoperatorFromRightOperator(Psinglet, PRsinglet) || !this->SuperoperatorFromRightOperator(Ptriplet, PRtriplet))
				return false;

			P = -1.0 * _operator->Rate1() * (PLsinglet * PRtriplet + PLtriplet * PRsinglet);

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationRandomFields)
		{
			// Adapted from Kattnig et al. (2016) DOI: 10.1088/1367-2630/18/6/063007

			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat E = arma::eye<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());	// Unity operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Rotate into the eigenbasis of the spin system
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from Sx
				if (!this->SuperoperatorFromOperators(Sx, Sx.t(), PB))
					return false;
				P += -1.0 * _operator->Rate1() * ((3 / 2) * E - PB);

				// Get the contribution from Sy
				if (!this->SuperoperatorFromOperators(Sy, Sy.t(), PB))
					return false;
				P += -1.0 * _operator->Rate2() * ((3 / 2) * E - PB);

				// Get the contribution from Sz
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB))
					return false;
				P += -1.0 * _operator->Rate3() * ((3 / 2) * E - PB);
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT1)
		{
			// Routine to calculate T1 relaxation (longitudinal)
			auto spins = _operator->Spins();

			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat S_plus;
			arma::sp_cx_mat S_minus;
			arma::sp_cx_mat PB; // Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL; // The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR; // The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sp(), (*i), S_plus);
				this->CreateOperator((*i)->Sm(), (*i), S_minus);

				// Rotate into the eigenbasis of the spin system
				S_plus = _rotationmatrix.t() * S_plus * _rotationmatrix;
				S_minus = _rotationmatrix.t() * S_minus * _rotationmatrix;

				// Get the contribution from T1 relaxation
				if (!this->SuperoperatorFromOperators(S_plus, S_minus.t(), PB) || !this->SuperoperatorFromLeftOperator(S_minus.t() * S_plus, PL) || !this->SuperoperatorFromRightOperator(S_minus.t() * S_plus, PR))
					return false;

				// \mathcal{L}_{T1}(\rho) = \frac{1}{T1} (2S_+ \rho S_- - \{S_-S_+, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::RelaxationT2)
		{
			// Routine to calculate T2 relaxation (transverse)
			auto spins = _operator->Spins();

			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::zeros<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions()); // Total operator
			arma::sp_cx_mat PB;																					// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;																					// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;																					// The right-hand operator (. * conj(P)*P)

			for (auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if (!this->Contains(*i))
					continue;

				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);

				// Rotate into the eigenbasis of the spin system
				Sx = _rotationmatrix.t() * Sx * _rotationmatrix;
				Sy = _rotationmatrix.t() * Sy * _rotationmatrix;
				Sz = _rotationmatrix.t() * Sz * _rotationmatrix;

				// Get the contribution from T2 relaxation
				if (!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				// \mathcal{L}_{T2}(\rho) = \frac{1}{T2} (2S_z \rho S_z - \{S_z^2, \rho\})
				P += (PB - (PL + PR) / 2.0) * _operator->Rate1();
			}

			// Set the resulting operator
			_out = P;
		}
		else if (_operator->Type() == OperatorType::Unspecified)
		{
			// Cannot construct operator if the type is not specified
			return false;
		}
		else
		{
			std::cout << "Cannot construct relaxation operator for Operator " << _operator->Name() << "!" << std::endl;
			return false;
		}

		return true;
	}
}