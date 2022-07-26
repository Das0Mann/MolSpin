/////////////////////////////////////////////////////////////////////////
// SpinSpace class (SpinAPI Module)
// ------------------
// This source file contains methods related to relaxation operators,
// which are described by Operator objects.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
namespace SpinAPI
{
	// -----------------------------------------------------
	// Relaxation operators
	// -----------------------------------------------------
	bool SpinSpace::RelaxationOperator(const operator_ptr& _operator, arma::cx_mat& _out) const
	{
		// Make sure that we have a valid operator object
		if(_operator == nullptr || !_operator->IsValid())
			return false;
		
		// The relaxation operators can only be constructed in the superspace
		if(!this->useSuperspace)
			return false;
		
		if(_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();
			arma::cx_mat Sx;
			arma::cx_mat Sy;
			arma::cx_mat Sz;
			arma::cx_mat P = arma::zeros<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());	// Total operator
			arma::cx_mat PB;	// Operator from both sides (P * . * conj(P))
			arma::cx_mat PL;	// The left-hand operator (conj(P)*P * .)
			arma::cx_mat PR;	// The right-hand operator (. * conj(P)*P)
			
			for(auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if(!this->Contains(*i))
					continue;
				
				// Create the spin operators
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sx()), (*i), Sx);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sy()), (*i), Sy);
				this->CreateOperator(arma::conv_to<arma::cx_mat>::from((*i)->Sz()), (*i), Sz);
				
				// Get the contribution from Sx
				if(!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR)/2.0) * _operator->Rate1();
				
				// Get the contribution from Sy
				if(!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR)/2.0) * _operator->Rate2();
				
				// Get the contribution from Sz
				if(!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR)/2.0) * _operator->Rate3();
			}
			
			// Set the resulting operator
			_out = P;
		}
		else if(_operator->Type() == OperatorType::Unspecified)
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
	bool SpinSpace::RelaxationOperator(const operator_ptr& _operator, arma::sp_cx_mat& _out) const
	{
		// Make sure that we have a valid operator object
		if(_operator == nullptr || !_operator->IsValid())
			return false;
		
		// The relaxation operators can only be constructed in the superspace
		if(!this->useSuperspace)
			return false;
		
		if(_operator->Type() == OperatorType::RelaxationLindblad)
		{
			auto spins = _operator->Spins();
			arma::sp_cx_mat Sx;
			arma::sp_cx_mat Sy;
			arma::sp_cx_mat Sz;
			arma::sp_cx_mat P = arma::sp_cx_mat(this->SpaceDimensions(), this->SpaceDimensions());	// Total operator
			arma::sp_cx_mat PB;	// Operator from both sides (P * . * conj(P))
			arma::sp_cx_mat PL;	// The left-hand operator (conj(P)*P * .)
			arma::sp_cx_mat PR;	// The right-hand operator (. * conj(P)*P)
			
			for(auto i = spins.cbegin(); i != spins.cend(); i++)
			{
				// Skip spins that are not part of the current SpinSpace
				if(!this->Contains(*i))
					continue;
				
				// Create the spin operators
				this->CreateOperator((*i)->Sx(), (*i), Sx);
				this->CreateOperator((*i)->Sy(), (*i), Sy);
				this->CreateOperator((*i)->Sz(), (*i), Sz);
				
				// Get the contribution from Sx
				if(!this->SuperoperatorFromOperators(Sx, Sx.t(), PB) || !this->SuperoperatorFromLeftOperator(Sx.t() * Sx, PL) || !this->SuperoperatorFromRightOperator(Sx.t() * Sx, PR))
					return false;
				P += (PB - (PL + PR)/2.0) * _operator->Rate1();
				
				// Get the contribution from Sy
				if(!this->SuperoperatorFromOperators(Sy, Sy.t(), PB) || !this->SuperoperatorFromLeftOperator(Sy.t() * Sy, PL) || !this->SuperoperatorFromRightOperator(Sy.t() * Sy, PR))
					return false;
				P += (PB - (PL + PR)/2.0) * _operator->Rate2();
				
				// Get the contribution from Sz
				if(!this->SuperoperatorFromOperators(Sz, Sz.t(), PB) || !this->SuperoperatorFromLeftOperator(Sz.t() * Sz, PL) || !this->SuperoperatorFromRightOperator(Sz.t() * Sz, PR))
					return false;
				P += (PB - (PL + PR)/2.0) * _operator->Rate3();
			}
			
			// Set the resulting operator
			_out = P;
		}
		else
		{
			std::cout << "Cannot construct relaxation operator for Operator " << _operator->Name() << "!" << std::endl;
			return false;
		}
		
		return true;
	}
}

