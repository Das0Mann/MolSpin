/////////////////////////////////////////////////////////////////////////
// SpinSpace class (SpinAPI Module)
// ------------------
// This source file contains methods for creating single-spin operators,
// and for converting operators or vectors to/from Superspace.
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
namespace SpinAPI
{
	// -----------------------------------------------------
	// Operator representations in the space
	// ---
	// Returns an operator (matrix) defined in the spin space.
	// Note that the matrix dimensions must fit NxN where N is the multiplicity of the spin.
	// -----------------------------------------------------
	// A block in the diagonal (dense)
	bool SpinSpace::CreateOperator(const arma::cx_mat& _operator, const spin_ptr& _spin, arma::cx_mat& _out) const
	{
		// Validate the input
		if(_spin == nullptr)
			return false;
		
		// Check that the given operator is valid
		if(_operator.n_rows != static_cast<unsigned int>(_spin->Multiplicity()) || _operator.n_cols != static_cast<unsigned int>(_spin->Multiplicity()))
			return false;
		
		// If the spin is outside the space, we just have an identity inside the space
		if(!this->Contains(_spin))
		{
			_out = arma::eye<arma::cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());
			return true;
		}
		
		// Prepare by getting the matrix in the space of the first spin
		arma::cx_mat result;
		auto i = this->spins.cbegin();
		
		if((*i) != _spin)
			result = arma::eye<arma::cx_mat>((*i)->Multiplicity(), (*i)->Multiplicity());
		else
			result = _operator;
		
		// Iterate over the remaining spins
		++i;	// We used the first one
		arma::cx_mat tmp;
		for( ; i != this->spins.cend(); i++)
		{
			if((*i) != _spin)
				tmp = arma::kron(result,arma::eye<arma::cx_mat>((*i)->Multiplicity(), (*i)->Multiplicity()));
			else
				tmp = arma::kron(result,_operator);
			
			result = tmp;
		}
		
		_out = result;
		return true;
	}
	
	// Create operator in the spin space from a single-spin operator (sparse version)
	bool SpinSpace::CreateOperator(const arma::sp_cx_mat& _operator, const spin_ptr& _spin, arma::sp_cx_mat& _out) const
	{
		// Validate the input
		if(_spin == nullptr)
			return false;
		
		// Check that the given operator is valid
		if(_operator.n_rows != static_cast<unsigned int>(_spin->Multiplicity()) || _operator.n_cols != static_cast<unsigned int>(_spin->Multiplicity()))
			return false;
		
		// If the spin is outside the space, we just have an identity inside the space
		if(!this->Contains(_spin))
		{
			_out = arma::speye<arma::sp_cx_mat>(this->SpaceDimensions(), this->SpaceDimensions());
			return true;
		}
		
		// Prepare by getting the matrix in the space of the first spin
		arma::sp_cx_mat result;
		auto i = this->spins.cbegin();
		
		if((*i) != _spin)
			result = arma::speye<arma::sp_cx_mat>((*i)->Multiplicity(), (*i)->Multiplicity());
		else
			result = _operator;
		
		// Iterate over the remaining spins
		++i;	// We used the first one
		arma::sp_cx_mat tmp;
		for( ; i != this->spins.cend(); i++)
		{
			if((*i) != _spin)
				tmp = arma::kron(result,arma::speye<arma::sp_cx_mat>((*i)->Multiplicity(), (*i)->Multiplicity()));
			else
				tmp = arma::kron(result,_operator);
			
			result = tmp;
		}
		
		_out = result;
		return true;
	}
	
	// Converts an operator to a superspace-vector
	bool SpinSpace::OperatorToSuperspace(const arma::cx_mat& _in, arma::cx_vec& _out) const	// TODO: Check matrix dimensions
	{
		_out = arma::vectorise(_in.st());
		return true;
	}
	
	// Sparse version
	// Implemented in terms of the dense version as it produces a "dense vector" anyway
	bool SpinSpace::OperatorToSuperspace(const arma::sp_cx_mat& _in, arma::cx_vec& _out) const
	{
		arma::cx_mat tmp(_in);
		return this->OperatorToSuperspace(tmp, _out);
	}

	// Converts a superspace vector to an operator in the Hilbert space
	bool SpinSpace::OperatorFromSuperspace(const arma::cx_vec& _in, arma::cx_mat& _out) const	// TODO: Check matrix dimensions
	{
		_out = arma::cx_mat(_in);
		_out.reshape(sqrt(_in.size()),sqrt(_in.size()));
		
		// Change by Luca Gerhards -- Needs to be done if matrix is not fully symmetric due to the filling of rows in c++
		_out = _out.st();
		
		return true;
	}
	
	// Sparse version
	// Implemented in terms of the dense version as it relies on a "dense vector" anyway
	bool SpinSpace::OperatorFromSuperspace(const arma::cx_vec& _in, arma::sp_cx_mat& _out) const
	{
		arma::cx_mat output;
		auto result = this->OperatorFromSuperspace(_in,output);
		_out = arma::conv_to<arma::sp_cx_mat>::from(output);
		return result;
	}
	
	// Defines the superoperator A(x) such that A(x) ~ L*x*R^T, where L and R are the left-side and right-side operators acting on an operator x, and ^T means transpose
	// Note: If you use an identity operator for either the left- or right-side operator, it is more efficient to use the other methods defined below.
	bool SpinSpace::SuperoperatorFromOperators(const arma::cx_mat& _leftside, const arma::cx_mat& _rightside, arma::cx_mat& _out) const	// TODO: Check matrix dimensions
	{
		_out = kron(_leftside, _rightside.st());
		return true;
	}
	
	// Sparse version
	bool SpinSpace::SuperoperatorFromOperators(const arma::sp_cx_mat& _leftside, const arma::sp_cx_mat& _rightside, arma::sp_cx_mat& _out) const
	{
		_out = kron(_leftside, _rightside.st());
		return true;
	}
		
	// Defines the superoperator A(x) such that A(x) ~ L*x*I, where L is the left-side operator acting on an operator x, and I is the identity
	bool SpinSpace::SuperoperatorFromLeftOperator(const arma::cx_mat& _leftside, arma::cx_mat& _out) const	// TODO: Implement more efficient algorithm, and check matrix dimensions
	{
		_out = kron(_leftside, arma::eye<arma::cx_mat>(arma::size(_leftside)));
		return true;
	}
	
	// Sparse version
	bool SpinSpace::SuperoperatorFromLeftOperator(const arma::sp_cx_mat& _leftside, arma::sp_cx_mat& _out) const
	{
		_out = kron(_leftside, arma::speye<arma::sp_cx_mat>(arma::size(_leftside)));
		return true;
	}
	
	// Defines the superoperator A(x) such that A(x) ~ I*x*R^T, where R is the right-side operator acting on an operator x, I is the identity, and ^T means transpose
	bool SpinSpace::SuperoperatorFromRightOperator(const arma::cx_mat& _rightside, arma::cx_mat& _out) const	// TODO: Implement more efficient algorithm, and check matrix dimensions
	{
		_out = kron(arma::eye<arma::cx_mat>(arma::size(_rightside)), _rightside.st());
		return true;
	}
	
	// Sparse version
	bool SpinSpace::SuperoperatorFromRightOperator(const arma::sp_cx_mat& _rightside, arma::sp_cx_mat& _out) const
	{
		_out = kron(arma::speye<arma::sp_cx_mat>(arma::size(_rightside)), _rightside.st());
		return true;
	}
	
	// -----------------------------------------------------
	// Methods for reordering the basis, i.e. change the
	// order in which the individual spins appear in the
	// Kronecker products.
	// -----------------------------------------------------
	// Reordering to the native order of the SpinSpace is done in terms of the more general method
	bool SpinSpace::ReorderBasis(arma::cx_vec& _out, const std::vector<spin_ptr>& _oldBasis) const
	{
		return this->ReorderBasis(_out, _oldBasis, this->spins);
	}
	bool SpinSpace::ReorderBasis(arma::cx_mat& _out, const std::vector<spin_ptr>& _oldBasis) const
	{
		return this->ReorderBasis(_out, _oldBasis, this->spins);
	}
	bool SpinSpace::ReorderBasis(arma::sp_cx_mat& _out, const std::vector<spin_ptr>& _oldBasis) const
	{
		return this->ReorderBasis(_out, _oldBasis, this->spins);
	}
	
	// Reordering by application of the reordering operator (vector)
	bool SpinSpace::ReorderBasis(arma::cx_vec& _out, const std::vector<spin_ptr>& _oldBasis, const std::vector<spin_ptr>& _newBasis) const
	{
		// Obtain the reordering operator
		arma::sp_cx_mat R;
		if(!this->ReorderingOperator(R, _oldBasis, _newBasis))
			return false;
		
		// Perform the reordering
		arma::cx_vec tmp = R * _out;
		_out = tmp;
		
		return true;
	}
	
	// Reordering by application of the reordering operator (dense matrix)
	bool SpinSpace::ReorderBasis(arma::cx_mat& _out, const std::vector<spin_ptr>& _oldBasis, const std::vector<spin_ptr>& _newBasis) const
	{
		// Obtain the reordering operator
		arma::sp_cx_mat R;
		if(!this->ReorderingOperator(R, _oldBasis, _newBasis))
			return false;
		
		// Perform the reordering
		arma::cx_mat tmp = R * _out * R.t();
		_out = tmp;
		
		return true;
	}
	
	// Reordering by application of the reordering operator (sparse matrix)
	bool SpinSpace::ReorderBasis(arma::sp_cx_mat& _out, const std::vector<spin_ptr>& _oldBasis, const std::vector<spin_ptr>& _newBasis) const
	{
		// Obtain the reordering operator
		arma::sp_cx_mat R;
		if(!this->ReorderingOperator(R, _oldBasis, _newBasis))
			return false;
		
		// Perform the reordering
		arma::sp_cx_mat tmp = R * _out * R.t();
		_out = tmp;
		
		return true;
	}
	
	// Calculation of the reordering operator
	bool SpinSpace::ReorderingOperator(arma::sp_cx_mat& _out, const std::vector<spin_ptr>& _oldBasis, const std::vector<spin_ptr>& _newBasis) const
	{
		// Check that the two spin collections contain the same number of spins (we are only reordering here!)
		if(_oldBasis.size() != _newBasis.size())
			return false;
		
		// The contained spins should also be identical
		for(auto i = _oldBasis.cbegin(); i != _oldBasis.cend(); i++)
			if(std::find(_newBasis.cbegin(), _newBasis.cend(), *i) == _newBasis.cend())
				return false;
		
		// Get the dimensions of the spin space
		unsigned int dimensions = 1;
		for(auto i = _oldBasis.cbegin(); i != _oldBasis.cend(); i++)
			dimensions *= static_cast<unsigned int>( (*i)->Multiplicity() );
		
		// Compute the reordering operator
		bool isFirst = true;
		std::vector<spin_ptr> tempBasis = _oldBasis;
		arma::sp_cx_mat R(dimensions, dimensions);
		arma::sp_cx_mat PrevId;
		arma::sp_cx_mat PostId;
		for(auto i = _newBasis.cbegin(); i != _newBasis.cend(); i++)
		{
			// Get an iterator to the current spin in the old basis, and its position in the old basis
			auto target = std::find(tempBasis.cbegin(), tempBasis.cend(), *i);
			auto targetPosition = target - tempBasis.cbegin();
			auto newPosition = i - _newBasis.cbegin();
			
			// Check whether both spins have the same position, as no action is required in that case
			while(targetPosition > newPosition)
			{
				// Get the dimensions of the spin subspace consiting of the spins before the target and its previous spin (that we are swapping with)
				unsigned int subsetDimensions = 1;
				for(auto j = tempBasis.cbegin(); j != target-1; j++)
					subsetDimensions *= static_cast<unsigned int>( (*j)->Multiplicity() );
				
				// Get the identity operator of the spin subspace that has been fixed so far
				PrevId = arma::eye<arma::sp_cx_mat>(subsetDimensions, subsetDimensions);
				
				// Also get the identity operator of the remaining subspace
				subsetDimensions = 1;
				for(auto j = target+1; j != tempBasis.cend(); j++)
					subsetDimensions *= static_cast<unsigned int>( (*j)->Multiplicity() );
				PostId = arma::eye<arma::sp_cx_mat>(subsetDimensions, subsetDimensions);
				
				// Produce swapping operator in the two-spin subspace of the spins to be swapped
				// NOTE: This operator is known as the perfect shuffle of the Kronecker product
				unsigned int p = static_cast<unsigned int>( (*(target-1))->Multiplicity() );
				unsigned int q = static_cast<unsigned int>( (*target)->Multiplicity() );
				unsigned int r = p*q;
				arma::sp_cx_mat S = arma::zeros<arma::sp_cx_mat>(r, r);
				unsigned int s = 0;
				for(unsigned int j = 0; j < q; j++)
				{
					for(unsigned int k = j; k < r; k += q)
						S(s++,k) = 1;
				}
				
				// Get the total swapping operator
				arma::sp_cx_mat P = kron(PrevId, kron(S, PostId));
				
				// Append the operator to the total reordering operator
				if(isFirst)
				{
					R = P;
					isFirst = false;
				}
				else
				{
					S = P * R;
					R = S;
				}
				
				// Swap neighbours in the old basis
				std::iter_swap(tempBasis.begin() + targetPosition, tempBasis.begin() + targetPosition - 1);
				
				// The target position is now closer to the beginning. Find is used again to make sure the iterator is valid after swap
				target = std::find(tempBasis.cbegin(), tempBasis.cend(), *i);
				targetPosition = target - tempBasis.cbegin();
			}
		}
		
		// If no reordering is necessary, return an identity
		if(isFirst)
			R = arma::eye<arma::sp_cx_mat>(dimensions, dimensions);
		
		// Set the result
		_out = R;
		
		return true;
	}
	
	// -----------------------------------------------------
	// Spherical tensors
	// -----------------------------------------------------
	

	//-----------------------------------------------------------------
	//-----------------RANK 1 TENSORS - BILINEAR-----------------------
	//-----------------------------------------------------------------

	// Rank 1 spherical tensors with m=0
	bool SpinSpace::Rk1SphericalTensorT0(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{

		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;
		arma::cx_mat S2x;
		arma::cx_mat S2y;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = -1.0 * (arma::cx_double(0.0, 1.0) / sqrt(2.0)) * ((S1x * S2y) - (S1y * S2x));	
		return true;
	}

	//Rank 1 spherical tensors with m=+1
	bool SpinSpace::Rk1SphericalTensorTp1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;
		arma::cx_mat S2x;
		arma::cx_mat S2y;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = -0.5 * (S1z * S2x - S1x * S2y + (arma::cx_double(0.0,1.0) * (S1z * S2y - S1y * S2z)));
		
		return true;
	}

	//Rank 1 spherical tensors with m=-1
	bool SpinSpace::Rk1SphericalTensorTm1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;
		arma::cx_mat S2x;
		arma::cx_mat S2y;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = -0.5 * (S1z * S2x - S1x * S2y + (arma::cx_double(0.0,1.0) * (S1z * S2y - S1y * S2z)));
		
		return true;
	}

	//-----------------------------------------------------------------
	//---------------------LINEAR RANK 2 TENSORS---------------------
	//-----------------------------------------------------------------

	// Rank 0 spherical tensor
	bool SpinSpace::LRk0TensorT0(const spin_ptr& _spin1,const arma::cx_vec& _field, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = (1.0/sqrt(3.0)) * (S1x * _field(0) + S1y * _field(1) + S1z * _field(2));
		return true;
	}

	// Rank 0 spherical tensor
	bool SpinSpace::LRk0TensorT0(const spin_ptr& _spin1,const arma::cx_vec& _field, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = (1.0/sqrt(3.0)) * (S1x * _field(0) + S1y * _field(1) + S1z * _field(2));
		return true;
	}

	// Rank 2 spherical tensor with m=0
	bool SpinSpace::LRk2SphericalTensorT0(const spin_ptr& _spin1,const arma::cx_vec& _field, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = (1.00/sqrt(6.0)) * ((3 * S1z * _field(2))-(S1x * _field(0) + S1y * _field(1) + S1z * _field(2)));

		//_out = (S1p*S2m + S1m*S2p + 2.0*S1z*S2z) * 0.40824829;	// Factor of 1/sqrt(6) = 0.40824829
		return true;
	}


	// Rank 2 spherical tensor with m=0, sparse version
	bool SpinSpace::LRk2SphericalTensorT0(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = (1.00/sqrt(6.0)) * ((3* S1z * _field(2))-(S1x * _field(0) + S1y * _field(1) + S1z * _field(2)));
		//_out = (S1p*S2m + S1m*S2p + 2.0*S1z*S2z) * 0.40824829;	// Factor of 1/sqrt(6) = 0.40824829
		return true;
	}
	
	// Rank 2 spherical tensor with m=+1
	bool SpinSpace::LRk2SphericalTensorTp1(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;

		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = -0.5 * (S1x * _field(2) + S1z * _field(0) + (arma::cx_double(0.0, 1.00) * (S1y * _field(2) + S1z * _field(1))));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=+1, sparse version
	bool SpinSpace::LRk2SphericalTensorTp1(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;

		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = -0.5 * (S1x * _field(2) + S1z * _field(0) + (arma::cx_double(0.0, 1.00) * (S1y * _field(2) + S1z * _field(1))));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=-1
	bool SpinSpace::LRk2SphericalTensorTm1(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;

		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * _field(2) + S1z * _field(0) - (arma::cx_double(0.0, 1.00) * (S1y * _field(2) + S1z * _field(1))));
		//_out = (S1m*S2z + S1z*S2m) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=-1, sparse version
	bool SpinSpace::LRk2SphericalTensorTm1(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;

		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * _field(2) + S1z * _field(0) - (arma::cx_double(0.0, 1.00) * (S1y * _field(2) + S1z * _field(1))));
		//_out = (S1m*S2z + S1z*S2m) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=+2
	bool SpinSpace::LRk2SphericalTensorTp2(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;

		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * _field(0) - S1y * _field(1) + (arma::cx_double(0.0, 1.0) * (S1x * _field(1) + S1y * _field(0))));
		//_out = S1p*S2p;
		return true;
	}
	
	// Rank 2 spherical tensor with m=+2, sparse version
	bool SpinSpace::LRk2SphericalTensorTp2(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;

		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * _field(0) - S1y * _field(1) + (arma::cx_double(0.0, 1.0) * (S1x * _field(1) + S1y * _field(0))));
		//_out = S1p*S2p;
		return true;
	}
	
	// Rank 2 spherical tensor with m=-2
	bool SpinSpace::LRk2SphericalTensorTm2(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;

		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * _field(0) - S1y * _field(1) - (arma::cx_double(0.0, 1.0) * (S1x * _field(1) + S1y * _field(0))));
		//_out = S1m*S2m;
		return true;
	}
	
	// Rank 2 spherical tensor with m=-2, sparse version
	bool SpinSpace::LRk2SphericalTensorTm2(const spin_ptr& _spin1, const arma::cx_vec& _field, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;

		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * _field(0) - S1y * _field(1) - (arma::cx_double(0.0, 1.0) * (S1x * _field(1) + S1y * _field(0))));
		//_out = S1m*S2m;
		return true;
	}

	//-----------------------------------------------------------------
	//---------------------BILINEAR RANK 2 TENSORS---------------------
	//-----------------------------------------------------------------
	
	// Rank 0 spherical tensor
	bool SpinSpace::BlRk0TensorT0(const spin_ptr& _spin1,const spin_ptr& _spin2 , arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S1y;
		arma::cx_mat S1z;
		arma::cx_mat S2x;
		arma::cx_mat S2y;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out =  (1.0/sqrt(3.0)) * (S1x * S2x + S1y * S2y + S1z * S2z );
		return true;
	}

	// Rank 0 spherical tensor, sparse version
	bool SpinSpace::BlRk0TensorT0(const spin_ptr& _spin1,const spin_ptr& _spin2 , arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S1z;
		arma::sp_cx_mat S2x;
		arma::sp_cx_mat S2y;
		arma::sp_cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out =  (1.0/sqrt(3.0)) * (S1x * S2x + S1y * S2y + S1z * S2z );
		return true;
	}

	// Rank 2 spherical tensor with m=0
	bool SpinSpace::BlRk2SphericalTensorT0(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S2x;
		arma::cx_mat S1y;
		arma::cx_mat S2y;
		arma::cx_mat S1z;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = (1.00/sqrt(6.00)) * (3.00 * S1z * S2z - (S1x * S2x + S1y * S2y + S1z * S2z));
		//_out = (S1p*S2m + S1m*S2p + 2.0*S1z*S2z) * 0.40824829;	// Factor of 1/sqrt(6) = 0.40824829
		return true;
	}


	// Rank 2 spherical tensor with m=0, sparse version
	bool SpinSpace::BlRk2SphericalTensorT0(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S2x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S2y;
		arma::sp_cx_mat S1z;
		arma::sp_cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}

		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = (1.00/sqrt(6.00)) * (3.00 * S1z * S2z - (S1x * S2x + S1y * S2y + S1z * S2z));
		//_out = (S1p*S2m + S1m*S2p + 2.0*S1z*S2z) * 0.40824829;	// Factor of 1/sqrt(6) = 0.40824829
		return true;
	}
	
	// Rank 2 spherical tensor with m=+1
	bool SpinSpace::BlRk2SphericalTensorTp1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S2x;
		arma::cx_mat S1y;
		arma::cx_mat S2y;
		arma::cx_mat S1z;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;
		// Create the spherical tensor
		_out = -0.5 * (S1x * S2z + S1z * S2x + (arma::cx_double(0.00, 1.00) * (S1y * S2z + S1z * S2y)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=+1, sparse version
	bool SpinSpace::BlRk2SphericalTensorTp1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S2x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S2y;
		arma::sp_cx_mat S1z;
		arma::sp_cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = -0.5 * (S1x * S2z + S1z * S2x + (arma::cx_double(0.00, 1.00) * (S1y * S2z + S1z * S2y)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=-1
	bool SpinSpace::BlRk2SphericalTensorTm1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S2x;
		arma::cx_mat S1y;
		arma::cx_mat S2y;
		arma::cx_mat S1z;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * S2z + S1z * S2x - (arma::cx_double(0.00, 1.00) * (S1y * S2z + S1z * S2y)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=-1, sparse version
	bool SpinSpace::BlRk2SphericalTensorTm1(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S2x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S2y;
		arma::sp_cx_mat S1z;
		arma::sp_cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * S2z + S1z * S2x - (arma::cx_double(0.00, 1.00) * (S1y * S2z + S1z * S2y)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=+2
	bool SpinSpace::BlRk2SphericalTensorTp2(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
	    arma::cx_mat S1x;
		arma::cx_mat S2x;
		arma::cx_mat S1y;
		arma::cx_mat S2y;
		arma::cx_mat S1z;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x *S2x - S1y *S2y + (arma::cx_double(0.0,1.0) * (S1x * S2y + S1y * S2x)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=+2, sparse version
	bool SpinSpace::BlRk2SphericalTensorTp2(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S2x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S2y;
		arma::sp_cx_mat S1z;
		arma::sp_cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x *S2x - S1y *S2y + (arma::cx_double(0.0,1.0) * (S1x * S2y + S1y * S2x)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=-2
	bool SpinSpace::BlRk2SphericalTensorTm2(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::cx_mat& _out) const
	{
		arma::cx_mat S1x;
		arma::cx_mat S2x;
		arma::cx_mat S1y;
		arma::cx_mat S2y;
		arma::cx_mat S1z;
		arma::cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x * S2x - S1y *S2y - (arma::cx_double(0.0,1.0) * (S1x * S2y + S1y * S2x)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
	
	// Rank 2 spherical tensor with m=-2, sparse version
	bool SpinSpace::BlRk2SphericalTensorTm2(const spin_ptr& _spin1, const spin_ptr& _spin2, arma::sp_cx_mat& _out) const
	{
		arma::sp_cx_mat S1x;
		arma::sp_cx_mat S2x;
		arma::sp_cx_mat S1y;
		arma::sp_cx_mat S2y;
		arma::sp_cx_mat S1z;
		arma::sp_cx_mat S2z;
		
		// We need all of these operators
		if(!this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sx() ), _spin1, S1x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sy() ), _spin1, S1y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin1->Sz() ), _spin1, S1z)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sx() ), _spin2, S2x)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sy() ), _spin2, S2y)
			|| !this->CreateOperator(arma::conv_to<arma::sp_cx_mat>::from( _spin2->Sz() ), _spin2, S2z))
		{
			return false;
		}
		
		S1x *= 2.0;
		S1y *= 2.0;
		S1z *= 2.0;
		S2x *= 2.0;
		S2y *= 2.0;
		S2z *= 2.0;

		// Create the spherical tensor
		_out = 0.5 * (S1x *S2x - S1y *S2y - (arma::cx_double(0.0,1.0) * (S1x * S2y + S1y * S2x)));
		//_out = (S1p*S2z + S1z*S2p) * 0.707106781;	// Factor of 1/sqrt(2) = 0.707106781
		return true;
	}
}

