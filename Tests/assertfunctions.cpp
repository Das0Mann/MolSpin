//////////////////////////////////////////////////////////////////////////////
// MolSpin Unit Testing Module
//
// Testing functions to check various conditions.
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
//////////////////////////////////////////////////////////////////////////////
// Compares two matrices to see whether they are equal
// Note: False for matrices with no elements!
// Real version
bool equal_matrices(const arma::mat &_m1, const arma::mat &_m2, double _tolerance = 1e-10)
{
	if (_m1.n_elem != _m2.n_elem || _m1.n_elem < 1)
		return false;

	auto diff = _m1 - _m2;

	return (arma::abs(diff).max() < _tolerance);
}

// Complex version in terms of the real version (real and imaginary parts taken separately)
bool equal_matrices(const arma::cx_mat &_m1, const arma::cx_mat &_m2, double _tolerance = 1e-10)
{
	return (equal_matrices(arma::real(_m1), arma::real(_m2), _tolerance) & equal_matrices(arma::imag(_m1), arma::imag(_m2), _tolerance));
}

// Sparse/dense variations
bool equal_matrices(const arma::sp_cx_mat &_m1, const arma::sp_cx_mat &_m2, double _tolerance = 1e-10)
{
	return equal_matrices(arma::conv_to<arma::cx_mat>::from(_m1), arma::conv_to<arma::cx_mat>::from(_m2), _tolerance);
}
bool equal_matrices(const arma::sp_cx_mat &_m1, const arma::cx_mat &_m2, double _tolerance = 1e-10)
{
	return equal_matrices(arma::conv_to<arma::cx_mat>::from(_m1), _m2, _tolerance);
}
bool equal_matrices(const arma::cx_mat &_m1, const arma::sp_cx_mat &_m2, double _tolerance = 1e-10)
{
	return equal_matrices(_m1, arma::conv_to<arma::cx_mat>::from(_m2), _tolerance);
}
//////////////////////////////////////////////////////////////////////////////
// Compares two floating point numbers
bool equal_double(double _d1, double _d2, double _tolerance = 1e-10)
{
	return (std::abs(_d1 - _d2) < _tolerance);
}
//////////////////////////////////////////////////////////////////////////////
// Compares two vectors
// Note: False for vectors with no elements!
bool equal_vec(arma::vec _v1, arma::vec _v2, double _tolerance = 1e-10)
{
	if (_v1.n_elem != _v2.n_elem || _v1.n_elem < 1)
		return false;

	return (norm(_v1 - _v2) < _tolerance);
}

// Complex vectors, implemented in terms of the real version
bool equal_vec(arma::cx_vec _v1, arma::cx_vec _v2, double _tolerance = 1e-10)
{
	return (equal_vec(arma::real(_v1), arma::real(_v2), _tolerance) & equal_vec(arma::imag(_v1), arma::imag(_v2), _tolerance));
}
//////////////////////////////////////////////////////////////////////////////
// Compares two tensors
bool equal_tensor(SpinAPI::Tensor _t1, SpinAPI::Tensor _t2, double _tolerance = 1e-10)
{
	bool isEqual = true;

	isEqual &= equal_double(_t1.Isotropic(), _t2.Isotropic(), _tolerance);
	isEqual &= equal_vec(_t1.Anisotropic(), _t2.Anisotropic(), _tolerance);
	isEqual &= equal_vec(_t1.Axis1(), _t2.Axis1(), _tolerance);
	isEqual &= equal_vec(_t1.Axis2(), _t2.Axis2(), _tolerance);
	isEqual &= equal_vec(_t1.Axis3(), _t2.Axis3(), _tolerance);
	isEqual &= equal_matrices(_t1.LabFrame(), _t2.LabFrame(), _tolerance);

	return isEqual;
}
//////////////////////////////////////////////////////////////////////////////
// Compares two strings with numbers, e.g. like "1.0 2.0 42"
bool equal_doublesfromstring(const std::string &_s1, const std::string &_s2, double _tolerance = 1e-10)
{
	std::vector<double> d1;
	std::vector<double> d2;
	d1.reserve(10);
	d2.reserve(10);
	std::string str = "";

	try
	{
		// Get numbers in the first string
		for (auto i = _s1.cbegin(); i != _s1.cend(); i++)
		{
			if (*i == '\t' || *i == '\n')
				continue;

			if (*i == ' ' && !str.empty())
			{
				d1.push_back(std::stod(str));
				str = "";
			}
			else
			{
				str += *i;
			}
		}

		// Get last element
		if (!str.empty())
		{
			d1.push_back(std::stod(str));
			str = "";
		}

		// Get numbers in the second string
		for (auto i = _s2.cbegin(); i != _s2.cend(); i++)
		{
			if (*i == '\t' || *i == '\n')
				continue;

			if (*i == ' ' && !str.empty())
			{
				d2.push_back(std::stod(str));
				str = "";
			}
			else
			{
				str += *i;
			}
		}

		// Get last element
		if (!str.empty())
		{
			d2.push_back(std::stod(str));
		}
	}
	catch (const std::exception &)
	{
		std::cout << "equal_doublesfromstring returns false due to an exception!";
		return false;
	}

	bool isEqual = true;

	// Check all the values
	auto id1 = d1.cbegin();
	auto id2 = d2.cbegin();
	while (id1 != d1.cend() && id2 != d2.cend())
	{
		isEqual &= equal_double(*id1, *id2);
		++id1;
		++id2;
	}

	// Check that both iterators has reached the end (i.e. same number of elements in both strings)
	isEqual &= (id1 == d1.cend() && id2 == d2.cend());

	return isEqual;
}
//////////////////////////////////////////////////////////////////////////////
