/////////////////////////////////////////////////////////////////////////
// Tensor class (SpinAPI Module)
// ------------------
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <sstream>
#include "Tensor.h"

namespace SpinAPI
{
	// -----------------------------------------------------
	// Spin Constructors and Destructor
	// -----------------------------------------------------
	Tensor::Tensor(double _isotropic) : isotropic(_isotropic), anisotropic(3, arma::fill::zeros), axis1({1, 0, 0}), axis2({0, 1, 0}), axis3({0, 0, 1}), mat1({0, 0, 0}), mat2({0, 0, 0}), mat3({0, 0, 0}), trajectory(),
										trjHasTime(false), trjHasIsotropic(false), trjHasAnisotropic(false), trjHasAxis1(false), trjHasAxis2(false), trjHasAxis3(false),
										trjHasMatXX(false), trjHasMatXY(false), trjHasMatXZ(false), trjHasMatYX(false), trjHasMatYY(false), trjHasMatYZ(false), trjHasMatZX(false), trjHasMatZY(false), trjHasMatZZ(false),
										trjTime(0), trjIsotropic(0), trjAnisotropicX(0), trjAnisotropicY(0), trjAnisotropicZ(0),
										trjAxis1X(0), trjAxis1Y(0), trjAxis1Z(0), trjAxis2X(0), trjAxis2Y(0), trjAxis2Z(0), trjAxis3X(0), trjAxis3Y(0), trjAxis3Z(0),
										trjMatXX(0), trjMatXY(0), trjMatXZ(0), trjMatYX(0), trjMatYY(0), trjMatYZ(0), trjMatZX(0), trjMatZY(0), trjMatZZ(0)
	{
	}

	Tensor::Tensor(double _isotropic, double _aniso1, double _aniso2, double _aniso3) : isotropic(_isotropic), anisotropic(3, arma::fill::zeros), axis1({1, 0, 0}), axis2({0, 1, 0}), axis3({0, 0, 1}), mat1({0, 0, 0}), mat2({0, 0, 0}), mat3({0, 0, 0}), trajectory(),
																						trjHasTime(false), trjHasIsotropic(false), trjHasAnisotropic(false), trjHasAxis1(false), trjHasAxis2(false), trjHasAxis3(false),
																						trjHasMatXX(false), trjHasMatXY(false), trjHasMatXZ(false), trjHasMatYX(false), trjHasMatYY(false), trjHasMatYZ(false), trjHasMatZX(false), trjHasMatZY(false), trjHasMatZZ(false),
																						trjTime(0), trjIsotropic(0), trjAnisotropicX(0), trjAnisotropicY(0), trjAnisotropicZ(0),
																						trjAxis1X(0), trjAxis1Y(0), trjAxis1Z(0), trjAxis2X(0), trjAxis2Y(0), trjAxis2Z(0), trjAxis3X(0), trjAxis3Y(0), trjAxis3Z(0),
																						trjMatXX(0), trjMatXY(0), trjMatXZ(0), trjMatYX(0), trjMatYY(0), trjMatYZ(0), trjMatZX(0), trjMatZY(0), trjMatZZ(0)
	{
		this->anisotropic(0) = _aniso1;
		this->anisotropic(1) = _aniso2;
		this->anisotropic(2) = _aniso3;
	}

	Tensor::Tensor(double _isotropic, const arma::vec &_anisotropic) : isotropic(_isotropic), anisotropic(_anisotropic), axis1({1, 0, 0}), axis2({0, 1, 0}), axis3({0, 0, 1}), mat1({0, 0, 0}), mat2({0, 0, 0}), mat3({0, 0, 0}), trajectory(),
																	   trjHasTime(false), trjHasIsotropic(false), trjHasAnisotropic(false), trjHasAxis1(false), trjHasAxis2(false), trjHasAxis3(false),
																	   trjHasMatXX(false), trjHasMatXY(false), trjHasMatXZ(false), trjHasMatYX(false), trjHasMatYY(false), trjHasMatYZ(false), trjHasMatZX(false), trjHasMatZY(false), trjHasMatZZ(false),
																	   trjTime(0), trjIsotropic(0), trjAnisotropicX(0), trjAnisotropicY(0), trjAnisotropicZ(0),
																	   trjAxis1X(0), trjAxis1Y(0), trjAxis1Z(0), trjAxis2X(0), trjAxis2Y(0), trjAxis2Z(0), trjAxis3X(0), trjAxis3Y(0), trjAxis3Z(0),
																	   trjMatXX(0), trjMatXY(0), trjMatXZ(0), trjMatYX(0), trjMatYY(0), trjMatYZ(0), trjMatZX(0), trjMatZY(0), trjMatZZ(0)
	{
	}

	Tensor::Tensor(double _isotropic, const arma::vec &_anisotropic, const arma::mat &_axes) : isotropic(_isotropic), anisotropic(_anisotropic), axis1(_axes.col(0)), axis2(_axes.col(1)), axis3(_axes.col(2)), mat1({0, 0, 0}), mat2({0, 0, 0}), mat3({0, 0, 0}), trajectory(),
																							   trjHasTime(false), trjHasIsotropic(false), trjHasAnisotropic(false), trjHasAxis1(false), trjHasAxis2(false), trjHasAxis3(false),
																							   trjHasMatXX(false), trjHasMatXY(false), trjHasMatXZ(false), trjHasMatYX(false), trjHasMatYY(false), trjHasMatYZ(false), trjHasMatZX(false), trjHasMatZY(false), trjHasMatZZ(false),
																							   trjTime(0), trjIsotropic(0), trjAnisotropicX(0), trjAnisotropicY(0), trjAnisotropicZ(0),
																							   trjAxis1X(0), trjAxis1Y(0), trjAxis1Z(0), trjAxis2X(0), trjAxis2Y(0), trjAxis2Z(0), trjAxis3X(0), trjAxis3Y(0), trjAxis3Z(0),
																							   trjMatXX(0), trjMatXY(0), trjMatXZ(0), trjMatYX(0), trjMatYY(0), trjMatYZ(0), trjMatZX(0), trjMatZY(0), trjMatZZ(0)
	{
		// Make sure axes are orthonormal
		this->OrthonormalizeAxes();
	}

	Tensor::Tensor(const arma::mat &_matrix) : isotropic(0.0), anisotropic(3, arma::fill::zeros), axis1({1, 0, 0}), axis2({0, 1, 0}), axis3({0, 0, 1}), mat1({0, 0, 0}), mat2({0, 0, 0}), mat3({0, 0, 0}), trajectory(),
											   trjHasTime(false), trjHasIsotropic(false), trjHasAnisotropic(false), trjHasAxis1(false), trjHasAxis2(false), trjHasAxis3(false),
											   trjHasMatXX(false), trjHasMatXY(false), trjHasMatXZ(false), trjHasMatYX(false), trjHasMatYY(false), trjHasMatYZ(false), trjHasMatZX(false), trjHasMatZY(false), trjHasMatZZ(false),
											   trjTime(0), trjIsotropic(0), trjAnisotropicX(0), trjAnisotropicY(0), trjAnisotropicZ(0),
											   trjAxis1X(0), trjAxis1Y(0), trjAxis1Z(0), trjAxis2X(0), trjAxis2Y(0), trjAxis2Z(0), trjAxis3X(0), trjAxis3Y(0), trjAxis3Z(0),
											   trjMatXX(0), trjMatXY(0), trjMatXZ(0), trjMatYX(0), trjMatYY(0), trjMatYZ(0), trjMatZX(0), trjMatZY(0), trjMatZZ(0)
	{
		// Diagonalize the matrix to obtain the principal axes and values
		this->DiagonalizeMatrix(_matrix);

		// And make sure the isotropic and anisotropic values are separated
		this->SeparateIsotropy();
	}

	Tensor::Tensor(const std::string &_tensor, const std::string &_path) : isotropic(0.0), anisotropic(3, arma::fill::zeros), axis1({1, 0, 0}), axis2({0, 1, 0}), axis3({0, 0, 1}), mat1({0, 0, 0}), mat2({0, 0, 0}), mat3({0, 0, 0}), trajectory(),
																		   trjHasTime(false), trjHasIsotropic(false), trjHasAnisotropic(false), trjHasAxis1(false), trjHasAxis2(false), trjHasAxis3(false),
																		   trjHasMatXX(false), trjHasMatXY(false), trjHasMatXZ(false), trjHasMatYX(false), trjHasMatYY(false), trjHasMatYZ(false), trjHasMatZX(false), trjHasMatZY(false), trjHasMatZZ(false),
																		   trjTime(0), trjIsotropic(0), trjAnisotropicX(0), trjAnisotropicY(0), trjAnisotropicZ(0),
																		   trjAxis1X(0), trjAxis1Y(0), trjAxis1Z(0), trjAxis2X(0), trjAxis2Y(0), trjAxis2Z(0), trjAxis3X(0), trjAxis3Y(0), trjAxis3Z(0),
																		   trjMatXX(0), trjMatXY(0), trjMatXZ(0), trjMatYX(0), trjMatYY(0), trjMatYZ(0), trjMatZX(0), trjMatZY(0), trjMatZZ(0)
	{
		if (!this->ParseTensor(_tensor, _path))
			std::cout << "Error: Failed to parse tensor!" << std::endl;
	}

	Tensor::Tensor(const Tensor &_tensor) : isotropic(_tensor.isotropic), anisotropic(_tensor.anisotropic), axis1(_tensor.axis1), axis2(_tensor.axis2), axis3(_tensor.axis3), mat1(_tensor.mat1), mat2(_tensor.mat2), mat3(_tensor.mat3), trajectory(_tensor.trajectory),
											trjHasTime(_tensor.trjHasTime), trjHasIsotropic(_tensor.trjHasIsotropic), trjHasAnisotropic(_tensor.trjHasAnisotropic), trjHasAxis1(_tensor.trjHasAxis1), trjHasAxis2(_tensor.trjHasAxis2), trjHasAxis3(_tensor.trjHasAxis3),
											trjHasMatXX(_tensor.trjHasMatXX), trjHasMatXY(_tensor.trjHasMatXY), trjHasMatXZ(_tensor.trjHasMatXZ), trjHasMatYX(_tensor.trjHasMatYX), trjHasMatYY(_tensor.trjHasMatYY), trjHasMatYZ(_tensor.trjHasMatYZ), trjHasMatZX(_tensor.trjHasMatZX), trjHasMatZY(_tensor.trjHasMatZY), trjHasMatZZ(_tensor.trjHasMatZZ),
											trjTime(_tensor.trjTime), trjIsotropic(_tensor.trjIsotropic), trjAnisotropicX(_tensor.trjAnisotropicX), trjAnisotropicY(_tensor.trjAnisotropicY), trjAnisotropicZ(_tensor.trjAnisotropicZ),
											trjAxis1X(_tensor.trjAxis1X), trjAxis1Y(_tensor.trjAxis1Y), trjAxis1Z(_tensor.trjAxis1Z), trjAxis2X(_tensor.trjAxis2X), trjAxis2Y(_tensor.trjAxis2Y), trjAxis2Z(_tensor.trjAxis2Z), trjAxis3X(_tensor.trjAxis3X), trjAxis3Y(_tensor.trjAxis3Y), trjAxis3Z(_tensor.trjAxis3Z),
											trjMatXX(_tensor.trjMatXX), trjMatXY(_tensor.trjMatXY), trjMatXZ(_tensor.trjMatXZ), trjMatYX(_tensor.trjMatYX), trjMatYY(_tensor.trjMatYY), trjMatYZ(_tensor.trjMatYZ), trjMatZX(_tensor.trjMatZX), trjMatZY(_tensor.trjMatZY), trjMatZZ(_tensor.trjMatZZ)
	{
	}

	Tensor::~Tensor()
	{
	}
	// -----------------------------------------------------
	// Operators
	// -----------------------------------------------------
	const Tensor &Tensor::operator=(const Tensor &_tensor)
	{
		this->isotropic = _tensor.isotropic;
		this->anisotropic = _tensor.anisotropic;
		this->axis1 = _tensor.axis1;
		this->axis2 = _tensor.axis2;
		this->axis3 = _tensor.axis3;
		this->trajectory = _tensor.trajectory;
		this->trjHasTime = _tensor.trjHasTime;
		this->trjHasIsotropic = _tensor.trjHasIsotropic;
		this->trjHasAnisotropic = _tensor.trjHasAnisotropic;
		this->trjHasAxis1 = _tensor.trjHasAxis1;
		this->trjHasAxis2 = _tensor.trjHasAxis2;
		this->trjHasAxis3 = _tensor.trjHasAxis3;

		this->trjHasMatXX = _tensor.trjHasMatXX;
		this->trjHasMatXY = _tensor.trjHasMatXY;
		this->trjHasMatXZ = _tensor.trjHasMatXZ;
		this->trjHasMatYX = _tensor.trjHasMatYX;
		this->trjHasMatYY = _tensor.trjHasMatYY;
		this->trjHasMatYZ = _tensor.trjHasMatYZ;
		this->trjHasMatZX = _tensor.trjHasMatZX;
		this->trjHasMatZY = _tensor.trjHasMatZY;
		this->trjHasMatZZ = _tensor.trjHasMatZZ;

		this->trjTime = _tensor.trjTime;
		this->trjIsotropic = _tensor.trjIsotropic;
		this->trjAnisotropicX = _tensor.trjAnisotropicX;
		this->trjAnisotropicY = _tensor.trjAnisotropicY;
		this->trjAnisotropicZ = _tensor.trjAnisotropicZ;
		this->trjAxis1X = _tensor.trjAxis1X;
		this->trjAxis1Y = _tensor.trjAxis1Y;
		this->trjAxis1Z = _tensor.trjAxis1Z;
		this->trjAxis2X = _tensor.trjAxis2X;
		this->trjAxis2Y = _tensor.trjAxis2Y;
		this->trjAxis2Z = _tensor.trjAxis2Z;
		this->trjAxis3X = _tensor.trjAxis3X;
		this->trjAxis3Y = _tensor.trjAxis3Y;
		this->trjAxis3Z = _tensor.trjAxis3Z;

		this->trjMatXX = _tensor.trjMatXX;
		this->trjMatXY = _tensor.trjMatXY;
		this->trjMatXZ = _tensor.trjMatXZ;
		this->trjMatYX = _tensor.trjMatYX;
		this->trjMatYY = _tensor.trjMatYY;
		this->trjMatYZ = _tensor.trjMatYZ;
		this->trjMatZX = _tensor.trjMatZX;
		this->trjMatZY = _tensor.trjMatZY;
		this->trjMatZZ = _tensor.trjMatZZ;

		return (*this);
	}

	// -----------------------------------------------------
	// Private methods
	// -----------------------------------------------------
	// Diagonalizes the matrix to obtain the principal axes and values
	void Tensor::DiagonalizeMatrix(const arma::mat &_matrix)
	{

		// Make sure that the matrix has the right dimensions
		if (_matrix.n_rows != 3 || _matrix.n_cols != 3)
		{
			std::cout << "Error: Cannot set Tensor by matrix of dimensions " << _matrix.n_rows << "x" << _matrix.n_cols << "! Must be 3x3! Ignoring matrix." << std::endl;
			return;
		}

		// Get a symmetric version of the matrix, i.e. the same matrix if is it symmetric
		arma::mat symmetrized_matrix = (_matrix + _matrix.t()) / 2.0;

		// Check whether the matrix was symmetric (i.e. _matrix - symmetrized_matrix == 0)
		if (abs(_matrix - symmetrized_matrix).max() > 1e-10)
		{
			std::cout << "Warning: Attempted to set Tensor from non-symmetric matrix." << std::endl;
			arma::cx_vec anisotropic_tmp;
			arma::cx_mat principalAxes_tmp = arma::eye<arma::cx_mat>(3, 3);
			;
			arma::eig_gen(anisotropic_tmp, principalAxes_tmp, _matrix);

			arma::mat principalAxes = arma::conv_to<arma::mat>::from(principalAxes_tmp);

			this->anisotropic = arma::conv_to<arma::vec>::from(anisotropic_tmp);
			this->axis1 = principalAxes.col(0);
			this->axis2 = principalAxes.col(1);
			this->axis3 = principalAxes.col(2);
		}
		else
		{
			// Do the diagonalization
			arma::mat principalAxes = arma::eye<arma::mat>(3, 3);
			arma::eig_sym(this->anisotropic, principalAxes, _matrix);
			this->axis1 = principalAxes.col(0);
			this->axis2 = principalAxes.col(1);
			this->axis3 = principalAxes.col(2);
		}
	}

	void Tensor::SeparateIsotropy()
	{
		// Put everything into the anisotropy vector
		this->anisotropic += this->isotropic;
		// Get the isotropic value
		this->isotropic = arma::sum(this->anisotropic) / 3.0;
		// And separate the isotropic value from the anisotropy
		this->anisotropic -= this->isotropic;
	}

	void Tensor::OrthonormalizeAxes()
	{
		// Make Axis1 and Axis2 orthogonal
		arma::vec tmp = dot(this->axis1, this->axis2) * this->axis1;
		this->axis2 -= tmp;

		// Make Axis3 orthogonal to Axis1 and Axis2
		tmp = dot(this->axis1, this->axis3) * this->axis1 + dot(this->axis2, this->axis3) * this->axis2;
		this->axis3 -= tmp;

		// Normalize the axes
		this->axis1 /= arma::norm(this->axis1);
		this->axis2 /= arma::norm(this->axis2);
		this->axis3 /= arma::norm(this->axis3);
	}

	double Tensor::DoubleFromString(std::string _str) const
	{
		std::istringstream stream(_str);
		double x;

		if (stream >> x)
			return x;

		return 0.0;
	}

	// Parse the given string to create a tensor
	// Used to parse tensor specifications from the input file
	bool Tensor::ParseTensor(const std::string &_str, const std::string &_path)
	{
		// Reset tensor
		this->isotropic = 0;
		this->anisotropic.zeros();
		bool hasTrajectory = false;

		// Split the input string at "+" to obtain a list of input elements
		std::vector<std::string> list;
		std::istringstream stream(_str);
		for (std::string str; std::getline(stream, str, '+'); list.push_back(str))
			;

		// Loop through the list
		for (auto i = list.cbegin(); i != list.cend(); i++)
		{
			auto j = i->find('(');
			if (j == std::string::npos)
			{
				// If just a number is given, parse that as the isotropic value
				this->isotropic += this->DoubleFromString(*i);
			}
			else
			{
				// Find the other paratheses
				auto k = i->find(')', j);
				if (k == std::string::npos)
				{
					// Parsing error; no matching parantheses
					return false;
				}

				// Obtain the keyword and value
				std::string keyword = i->substr(0, j);
				std::string value = i->substr(j + 1, k - j - 1);

				// Get the data
				if (keyword.compare("isotropic") == 0)
				{
					this->isotropic += this->DoubleFromString(value);
				}
				else if (keyword.compare("anisotropic") == 0)
				{
					arma::vec aniso(value);

					// Make sure that the anisotropic part has exactly 3 elements
					if (aniso.n_elem != 3)
					{
						std::cout << "Error: A tensor must have 3 values for the anisotropic part! " << aniso.n_elem << " values was given." << std::endl;
						return false;
					}

					this->anisotropic += aniso;
				}
				else if (keyword.compare("axis1") == 0)
				{
					this->axis1 = arma::vec(value);
				}
				else if (keyword.compare("axis2") == 0)
				{
					this->axis2 = arma::vec(value);
				}
				else if (keyword.compare("axis3") == 0)
				{
					this->axis3 = arma::vec(value);
				}
				else if (keyword.compare("matrix") == 0)
				{
					arma::mat tmp(value);
					arma::vec aniso(this->anisotropic); // Save the current value for the anisotropic part (can also contain isotropic contributions at this point)
					this->DiagonalizeMatrix(tmp);		// Sets anisotropic part and axes
					this->anisotropic += aniso;			// Add whatever anisotropies that were specified previously
				}
				else if (keyword.compare("changebasis") == 0)
				{
					arma::mat tmp(value);					   // Get the new basis
					arma::mat tmpLF = this->LabFrame();		   // Get the matrix representation in the standard basis
					arma::mat changed = tmp.t() * tmpLF * tmp; // Similarity transformation, note that matrix inversion can fail (if singular)
					this->DiagonalizeMatrix(changed);		   // Calculate a new tensor based on a similarity transformation
				}
				else if (keyword.compare("trajectory") == 0)
				{
					// Only a single trajectory can be specified
					if (hasTrajectory)
					{
						std::cout << "Error: Only a single trajectory can be assigned to a tensor!" << std::endl;
						return false;
					}

					hasTrajectory = true;
					this->LoadTrajectory(value, _path);
				}
				else
				{
					// A keyword was not recognized
					std::cout << "Error: Keyword \"" << keyword << "\" was not recognized when trying to parse tensor!" << std::endl;
					return false;
				}
			}
		}

		// Put the tensor data members into the correct format
		this->SeparateIsotropy();

		// Make sure the axes are orthonormal
		this->OrthonormalizeAxes();

		return true;
	}
	// -----------------------------------------------------
	// Methods to create ActionTarget objects
	// -----------------------------------------------------
	// Create ActionVectors
	std::vector<RunSection::NamedActionVector> Tensor::CreateActionVectors(const std::string &_system)
	{
		std::vector<RunSection::NamedActionVector> vectors;

		// Create ActionVectors
		RunSection::ActionVector anisotropicPart = RunSection::ActionVector(this->anisotropic, nullptr, true);
		RunSection::ActionVector axis1av = RunSection::ActionVector(this->axis1, nullptr, true);
		RunSection::ActionVector axis2av = RunSection::ActionVector(this->axis2, nullptr, true);
		RunSection::ActionVector axis3av = RunSection::ActionVector(this->axis3, nullptr, true);
		vectors.push_back(RunSection::NamedActionVector(_system + ".tensor.anisotropic", anisotropicPart));
		vectors.push_back(RunSection::NamedActionVector(_system + ".tensor.axis1", axis1av));
		vectors.push_back(RunSection::NamedActionVector(_system + ".tensor.axis2", axis2av));
		vectors.push_back(RunSection::NamedActionVector(_system + ".tensor.axis3", axis3av));

		return vectors;
	}

	// Create ActionScalars
	std::vector<RunSection::NamedActionScalar> Tensor::CreateActionScalars(const std::string &_system)
	{
		// Whether we have a trajectory; if so the ActionTargets should be readonly
		bool hasTrj = this->TrajectoryLength() != 0;

		std::vector<RunSection::NamedActionScalar> scalars;

		// Create ActionScalars
		RunSection::ActionScalar isotropicPart = RunSection::ActionScalar(this->isotropic, &CheckActionScalarTensorIsotropicPart, hasTrj);
		scalars.push_back(RunSection::NamedActionScalar(_system + ".tensor.isotropic", isotropicPart));

		return scalars;
	}

	// Method that calls the methods to generate ActionVectors and ActionScalars and inserts them into the given collections
	// Note that there are currently no ActionVectors as they cannot be readily added to this class
	void Tensor::GetActionTargets(std::vector<RunSection::NamedActionScalar> &_scalars, std::vector<RunSection::NamedActionVector> &_vectors, const std::string &_system)
	{
		// Get ActionTargets from private methods
		auto scalars = this->CreateActionScalars(_system);
		auto vectors = this->CreateActionVectors(_system);

		// Insert them
		_scalars.insert(_scalars.end(), scalars.begin(), scalars.end());
		_vectors.insert(_vectors.end(), vectors.begin(), vectors.end());
	}
	// -----------------------------------------------------
	// Public methods
	// -----------------------------------------------------
	// Return the isotropic part of the tensor
	double Tensor::Isotropic() const
	{
		return this->isotropic;
	}

	// Return the three anisotropic values
	arma::vec Tensor::Anisotropic() const
	{
		return this->anisotropic;
	}

	// Return the first principal axis (eigenvector)
	arma::vec Tensor::Axis1() const
	{
		return this->axis1;
	}

	// Return the second principal axis (eigenvector)
	arma::vec Tensor::Axis2() const
	{
		return this->axis2;
	}

	// Return the third principal axis (eigenvector)
	arma::vec Tensor::Axis3() const
	{
		return this->axis3;
	}

	// Return a matrix where each column is a principal axis (eigenvector)
	arma::mat Tensor::Axes() const
	{
		arma::mat axes = arma::eye<arma::mat>(3, 3);
		axes.col(0) = this->axis1;
		axes.col(1) = this->axis2;
		axes.col(2) = this->axis3;

		return axes;
	}

	// Returns the matrix representation of the tensor in the laboratory frame
	arma::mat Tensor::LabFrame() const
	{
		arma::mat axes = arma::eye<arma::mat>(3, 3);
		axes.col(0) = this->axis1;
		axes.col(1) = this->axis2;
		axes.col(2) = this->axis3;

		// TODO: Consider whether arma::inv(axes) should be used instead of axes.t() if "axes" are not orthogonal
		return (axes * arma::diagmat(this->anisotropic) * axes.t()) + arma::eye(size(axes)) * this->isotropic;
	}

	// Return the length of the trajectory (0 if no trajectory is assigned)
	unsigned int Tensor::TrajectoryLength() const
	{
		return this->trajectory.Length();
	}

	// Set the current state of the tensor to a specific step in the trajectory (row in the trajectory file)
	bool Tensor::SetTrajectoryStep(unsigned int _step)
	{
		// Make sure that we have a trajectory
		if (this->trajectory.Length() < 1)
			return false;

		// Make sure that the requested step exists
		if (_step >= this->trajectory.Length())
			return false;

		// Get the isotropic part
		if (this->trjHasIsotropic)
			this->isotropic = this->trajectory.Get(_step, this->trjIsotropic);

		// Get the anisotropic part
		if (this->trjHasAnisotropic)
		{
			this->anisotropic(0) = this->trajectory.Get(_step, this->trjAnisotropicX);
			this->anisotropic(1) = this->trajectory.Get(_step, this->trjAnisotropicY);
			this->anisotropic(2) = this->trajectory.Get(_step, this->trjAnisotropicZ);
			this->SeparateIsotropy();
		}

		// Get the axes
		if (this->trjHasAxis1)
		{
			this->axis1(0) = this->trajectory.Get(_step, this->trjAxis1X);
			this->axis1(1) = this->trajectory.Get(_step, this->trjAxis1Y);
			this->axis1(2) = this->trajectory.Get(_step, this->trjAxis1Z);
		}
		if (this->trjHasAxis2)
		{
			this->axis2(0) = this->trajectory.Get(_step, this->trjAxis2X);
			this->axis2(1) = this->trajectory.Get(_step, this->trjAxis2Y);
			this->axis2(2) = this->trajectory.Get(_step, this->trjAxis2Z);
		}
		if (this->trjHasAxis3)
		{
			this->axis3(0) = this->trajectory.Get(_step, this->trjAxis3X);
			this->axis3(1) = this->trajectory.Get(_step, this->trjAxis3Y);
			this->axis3(2) = this->trajectory.Get(_step, this->trjAxis3Z);
		}

		// Make sure that the axes are still orthonormal
		if (this->trjHasAxis1 || this->trjHasAxis2 || this->trjHasAxis3)
			this->OrthonormalizeAxes();

		// Get matrix elements
		if (this->trjHasMatXX)
		{
			this->mat1(0) = this->trajectory.Get(_step, this->trjMatXX);
		}

		if (this->trjHasMatXY)
		{
			this->mat1(1) = this->trajectory.Get(_step, this->trjMatXY);
		}
		if (this->trjHasMatXZ)
		{
			this->mat1(2) = this->trajectory.Get(_step, this->trjMatXZ);
		}

		if (this->trjHasMatYX)
		{
			this->mat2(0) = this->trajectory.Get(_step, this->trjMatYX);
		}
		if (this->trjHasMatYY)
		{
			this->mat2(1) = this->trajectory.Get(_step, this->trjMatYY);
		}
		if (this->trjHasMatYZ)
		{
			this->mat2(2) = this->trajectory.Get(_step, this->trjMatYZ);
		}

		if (this->trjHasMatZX)
		{
			this->mat3(0) = this->trajectory.Get(_step, this->trjMatZX);
		}
		if (this->trjHasMatZY)
		{
			this->mat3(1) = this->trajectory.Get(_step, this->trjMatZY);
		}
		if (this->trjHasMatZZ)
		{
			this->mat3(2) = this->trajectory.Get(_step, this->trjMatZZ);
		}

		// Make matrix diagonal to further use it in the code
		if (this->trjHasMatXX || this->trjHasMatXY || this->trjHasMatXZ || this->trjHasMatYX || this->trjHasMatYY || this->trjHasMatYZ || this->trjHasMatZX || this->trjHasMatZY || this->trjHasMatZZ)
		{

			arma::mat tmp_trj_mat = {{this->mat1(0), this->mat1(1), this->mat1(2)},
									 {this->mat2(0), this->mat2(1), this->mat2(2)},
									 {this->mat3(0), this->mat3(1), this->mat3(2)}};
			this->DiagonalizeMatrix(tmp_trj_mat); // Sets anisotropic part and axes
		}

		return true;
	}

	// Set the current state of the tensor to a specific time defined in the trajectory (requires a "time" column in the trajectory)
	// Uses linear interpolation based on the time
	bool Tensor::SetTime(double _time)
	{
		// The time must be specified in the trajectory in order to use this method
		if (!this->trjHasTime)
			return false;

		// Get the row corresponding to the time right after (or at) the requested time
		unsigned int row = 0;
		double timestepAbove = 0.0;
		if (!this->trajectory.FirstRowEqGreaterThan(_time, this->trjTime, row, timestepAbove))
		{
			// If no such time was found, i.e. we are past the last element of the trajectory, use the last step of the trajectory
			this->SetTrajectoryStep(this->trajectory.Length() - 1);
		}
		else if (row == 0)
		{
			// If the specified time is before the first element of the trajectory, use the first step of the trajectory
			this->SetTrajectoryStep(0);
		}
		else
		{
			double timestepBelow = this->trajectory.Get(row - 1, trjTime);			  // Get time at previous step
			double l = 1 - (timestepAbove - _time) / (timestepAbove - timestepBelow); // Get linear interpolation factor, i.e. a number between 0 and 1, where 0 = previous time step and 1 = next time step in the trajectory

			// Get the isotropic part
			if (this->trjHasIsotropic)
				this->isotropic = this->trajectory.Get(row, this->trjIsotropic) * l + this->trajectory.Get(row - 1, this->trjIsotropic) * (1 - l);

			// Get the anisotropic part
			if (this->trjHasAnisotropic)
			{
				this->anisotropic(0) = this->trajectory.Get(row, this->trjAnisotropicX) * l + this->trajectory.Get(row - 1, this->trjAnisotropicX) * (1 - l);
				this->anisotropic(1) = this->trajectory.Get(row, this->trjAnisotropicY) * l + this->trajectory.Get(row - 1, this->trjAnisotropicY) * (1 - l);
				this->anisotropic(2) = this->trajectory.Get(row, this->trjAnisotropicZ) * l + this->trajectory.Get(row - 1, this->trjAnisotropicZ) * (1 - l);
				this->SeparateIsotropy();
			}

			// Get the axes
			if (this->trjHasAxis1)
			{
				this->axis1(0) = this->trajectory.Get(row, this->trjAxis1X) * l + this->trajectory.Get(row - 1, this->trjAxis1X) * (1 - l);
				this->axis1(1) = this->trajectory.Get(row, this->trjAxis1Y) * l + this->trajectory.Get(row - 1, this->trjAxis1Y) * (1 - l);
				this->axis1(2) = this->trajectory.Get(row, this->trjAxis1Z) * l + this->trajectory.Get(row - 1, this->trjAxis1Z) * (1 - l);
			}
			if (this->trjHasAxis2)
			{
				this->axis2(0) = this->trajectory.Get(row, this->trjAxis2X) * l + this->trajectory.Get(row - 1, this->trjAxis2X) * (1 - l);
				this->axis2(1) = this->trajectory.Get(row, this->trjAxis2Y) * l + this->trajectory.Get(row - 1, this->trjAxis2Y) * (1 - l);
				this->axis2(2) = this->trajectory.Get(row, this->trjAxis2Z) * l + this->trajectory.Get(row - 1, this->trjAxis2Z) * (1 - l);
			}
			if (this->trjHasAxis3)
			{
				this->axis3(0) = this->trajectory.Get(row, this->trjAxis3X) * l + this->trajectory.Get(row - 1, this->trjAxis3X) * (1 - l);
				this->axis3(1) = this->trajectory.Get(row, this->trjAxis3Y) * l + this->trajectory.Get(row - 1, this->trjAxis3Y) * (1 - l);
				this->axis3(2) = this->trajectory.Get(row, this->trjAxis3Z) * l + this->trajectory.Get(row - 1, this->trjAxis3Z) * (1 - l);
			}

			// Make sure that the axes are still orthonormal
			if (this->trjHasAxis1 || this->trjHasAxis2 || this->trjHasAxis3)
				this->OrthonormalizeAxes();

			// get matrix elements
			if (this->trjHasMatXX)
			{
				this->mat1(0) = this->trajectory.Get(row, this->trjMatXX) * l + this->trajectory.Get(row - 1, this->trjMatXX) * (1 - l);
			}

			if (this->trjHasMatXY)
			{
				this->mat1(1) = this->trajectory.Get(row, this->trjMatXY) * l + this->trajectory.Get(row - 1, this->trjMatXY) * (1 - l);
			}
			if (this->trjHasMatXZ)
			{
				this->mat1(2) = this->trajectory.Get(row, this->trjMatXZ) * l + this->trajectory.Get(row - 1, this->trjMatXZ) * (1 - l);
			}

			if (this->trjHasMatYX)
			{
				this->mat2(0) = this->trajectory.Get(row, this->trjMatYX) * l + this->trajectory.Get(row - 1, this->trjMatYX) * (1 - l);
			}
			if (this->trjHasMatYY)
			{
				this->mat2(1) = this->trajectory.Get(row, this->trjMatYY) * l + this->trajectory.Get(row - 1, this->trjMatYY) * (1 - l);
			}
			if (this->trjHasMatYZ)
			{
				this->mat2(2) = this->trajectory.Get(row, this->trjMatYZ) * l + this->trajectory.Get(row - 1, this->trjMatYZ) * (1 - l);
			}

			if (this->trjHasMatZX)
			{
				this->mat3(0) = this->trajectory.Get(row, this->trjMatZX) * l + this->trajectory.Get(row - 1, this->trjMatZX) * (1 - l);
			}
			if (this->trjHasMatZY)
			{
				this->mat3(1) = this->trajectory.Get(row, this->trjMatZY) * l + this->trajectory.Get(row - 1, this->trjMatZY) * (1 - l);
			}
			if (this->trjHasMatZZ)
			{
				this->mat3(2) = this->trajectory.Get(row, this->trjMatZZ) * l + this->trajectory.Get(row - 1, this->trjMatZZ) * (1 - l);
			}

			// Make matrix diagonal to further use it in the code
			if (this->trjHasMatXX || this->trjHasMatXY || this->trjHasMatXZ || this->trjHasMatYX || this->trjHasMatYY || this->trjHasMatYZ || this->trjHasMatZX || this->trjHasMatZY || this->trjHasMatZZ)
			{

				arma::mat tmp_trj_mat = {{this->mat1(0), this->mat1(1), this->mat1(2)},
										 {this->mat2(0), this->mat2(1), this->mat2(2)},
										 {this->mat3(0), this->mat3(1), this->mat3(2)}};
				this->DiagonalizeMatrix(tmp_trj_mat); // Sets anisotropic part and axes
			}
		}

		return true;
	}

	// Used to set the necessary variables for time-dependent tensors specified by a tensortype
	void Tensor::SetTensor(arma::mat &m)
	{
		this->DiagonalizeMatrix(m);							  // diagonalise the matrix m - resets the anisotropic part
		this->isotropic = arma::sum(this->anisotropic) / 3.0; // separate the anisotropic and isotropic parts
		this->anisotropic -= this->isotropic;
	}

	// Returns lab frame representation of the tensor
	arma::mat Tensor::GetTensor()
	{
		arma::mat tensor = this->LabFrame();
		return tensor;
	}

	// Loads a trajectory and checks for headers used by the Tensor class
	bool Tensor::LoadTrajectory(const std::string &_filename, const std::string &_path, bool _overwrite)
	{
		bool loadedNewTrajectory = false;

		// Attempt to load new trajectory, but if a trajectory was already loaded check whether we can overwrite that
		if (this->trajectory.Length() < 1 || _overwrite)
			loadedNewTrajectory = this->trajectory.Load(_path + _filename);

		// If a new trajectory was loaded
		if (loadedNewTrajectory)
		{
			// Get time column
			this->trjHasTime = this->trajectory.HasColumn("time", trjTime);

			// Get isotropic value
			this->trjHasIsotropic = this->trajectory.HasColumn("isotropic", trjIsotropic);

			// Get anisotropic values
			this->trjHasAnisotropic = this->trajectory.HasColumn("anisotropic.x", trjAnisotropicX);
			this->trjHasAnisotropic &= this->trajectory.HasColumn("anisotropic.y", trjAnisotropicY);
			this->trjHasAnisotropic &= this->trajectory.HasColumn("anisotropic.z", trjAnisotropicZ);

			// Get axes
			this->trjHasAxis1 = this->trajectory.HasColumn("axis1.x", trjAxis1X);
			this->trjHasAxis1 &= this->trajectory.HasColumn("axis1.y", trjAxis1Y);
			this->trjHasAxis1 &= this->trajectory.HasColumn("axis1.z", trjAxis1Z);
			this->trjHasAxis2 = this->trajectory.HasColumn("axis2.x", trjAxis2X);
			this->trjHasAxis2 &= this->trajectory.HasColumn("axis2.y", trjAxis2Y);
			this->trjHasAxis2 &= this->trajectory.HasColumn("axis2.z", trjAxis2Z);
			this->trjHasAxis3 = this->trajectory.HasColumn("axis3.x", trjAxis3X);
			this->trjHasAxis3 &= this->trajectory.HasColumn("axis3.y", trjAxis3Y);
			this->trjHasAxis3 &= this->trajectory.HasColumn("axis3.z", trjAxis3Z);

			// matrix elements
			this->trjHasMatXX = this->trajectory.HasColumn("mat.xx", trjMatXX);
			this->trjHasMatXY = this->trajectory.HasColumn("mat.xy", trjMatXY);
			this->trjHasMatXZ = this->trajectory.HasColumn("mat.xz", trjMatXZ);
			this->trjHasMatYX = this->trajectory.HasColumn("mat.yx", trjMatYX);
			this->trjHasMatYY = this->trajectory.HasColumn("mat.yy", trjMatYY);
			this->trjHasMatYZ = this->trajectory.HasColumn("mat.yz", trjMatYZ);
			this->trjHasMatZX = this->trajectory.HasColumn("mat.zx", trjMatZX);
			this->trjHasMatZY = this->trajectory.HasColumn("mat.zy", trjMatZY);
			this->trjHasMatZZ = this->trajectory.HasColumn("mat.zz", trjMatZZ);

			// Check whether any useful data was contained in the trajectory (note that time is irrelevant here, as it is just used for indexing)
			if (!(this->trjHasIsotropic | this->trjHasAnisotropic | this->trjHasAxis1 | this->trjHasAxis2 | this->trjHasAxis3 | this->trjHasMatXX | this->trjHasMatXY | this->trjHasMatXZ | this->trjHasMatYX | this->trjHasMatYY | this->trjHasMatYZ | this->trjHasMatZX | this->trjHasMatZY | this->trjHasMatZZ))
			{
				// Unload the trajectory if it contained no useful information
				this->trajectory.Clear();
				loadedNewTrajectory = false;
			}
		}

		return loadedNewTrajectory;
	}
	// -----------------------------------------------------
	// Spherical tensor representation components
	// These are the components that should multiply the
	// various spherical tensor operators.
	// Note: The tensor is traceless (A_xx + A_yy + A_zz = 0) as the isotropic part is not included
	// Note: The tensor is symmetric, A_ij = A_ji
	// -----------------------------------------------------
	// T0 component: sqrt(6)/2 * A_zz
	arma::cx_double Tensor::SphericalT0() const
	{
		// Get the lab-frame A_zz component
		double A_zz = this->axis3(0) * this->axis3(0) * anisotropic(0);
		A_zz += this->axis3(1) * this->axis3(1) * anisotropic(1);
		A_zz += this->axis3(2) * this->axis3(2) * anisotropic(2);

		return arma::cx_double(1.224744871 * A_zz, 0.0); // Factor of sqrt(6)/2 = 1.224744871
	}

	// T+1 component: (A_xz + i*A_yz) / sqrt(2)
	arma::cx_double Tensor::SphericalTp1() const
	{
		// Get the lab-frame A_xz component
		double A_xz = this->axis1(0) * this->axis3(0) * anisotropic(0);
		A_xz += this->axis1(1) * this->axis3(1) * anisotropic(1);
		A_xz += this->axis1(2) * this->axis3(2) * anisotropic(2);

		// Get the lab-frame A_yz component
		double A_yz = this->axis2(0) * this->axis3(0) * anisotropic(0);
		A_yz += this->axis2(1) * this->axis3(1) * anisotropic(1);
		A_yz += this->axis2(2) * this->axis3(2) * anisotropic(2);

		return arma::cx_double(A_xz, A_yz) * 0.707106781; // Factor of 1/sqrt(2) = 0.707106781
	}

	// T-1 component: (A_xz - i*A_yz) / sqrt(2)
	arma::cx_double Tensor::SphericalTm1() const
	{
		// Get the lab-frame A_xz component
		double A_xz = this->axis1(0) * this->axis3(0) * anisotropic(0);
		A_xz += this->axis1(1) * this->axis3(1) * anisotropic(1);
		A_xz += this->axis1(2) * this->axis3(2) * anisotropic(2);

		// Get the lab-frame A_yz component
		double A_yz = this->axis2(0) * this->axis3(0) * anisotropic(0);
		A_yz += this->axis2(1) * this->axis3(1) * anisotropic(1);
		A_yz += this->axis2(2) * this->axis3(2) * anisotropic(2);

		return arma::cx_double(A_xz, -A_yz) * 0.707106781; // Factor of 1/sqrt(2) = 0.707106781
	}

	// T+2 component: (A_xx - A_yy) / 4 + (A_xy + A_yx) / 4i
	// Due to symmetry: (A_xx - A_yy) / 4 - i*A_xy / 2
	arma::cx_double Tensor::SphericalTp2() const
	{
		// Get the lab-frame A_xy component (equal to A_yx)
		double A_xy = this->axis1(0) * this->axis2(0) * anisotropic(0);
		A_xy += this->axis1(1) * this->axis2(1) * anisotropic(1);
		A_xy += this->axis1(2) * this->axis2(2) * anisotropic(2);

		// Get the lab-frame A_xx component
		double A_xx = this->axis1(0) * this->axis1(0) * anisotropic(0);
		A_xx += this->axis1(1) * this->axis1(1) * anisotropic(1);
		A_xx += this->axis1(2) * this->axis1(2) * anisotropic(2);

		// Get the lab-frame A_yy component
		double A_yy = this->axis2(0) * this->axis2(0) * anisotropic(0);
		A_yy += this->axis2(1) * this->axis2(1) * anisotropic(1);
		A_yy += this->axis2(2) * this->axis2(2) * anisotropic(2);

		return arma::cx_double((A_xx - A_yy) / 4.0, -A_xy / 2.0);
	}

	// T-2 component: (A_xx - A_yy) / 4 - (A_xy + A_yx) / 4i
	// Due to symmetry: (A_xx - A_yy) / 4 + i*A_xy / 2
	arma::cx_double Tensor::SphericalTm2() const
	{
		// Get the lab-frame A_xy component (equal to A_yx)
		double A_xy = this->axis1(0) * this->axis2(0) * anisotropic(0);
		A_xy += this->axis1(1) * this->axis2(1) * anisotropic(1);
		A_xy += this->axis1(2) * this->axis2(2) * anisotropic(2);

		// Get the lab-frame A_xx component
		double A_xx = this->axis1(0) * this->axis1(0) * anisotropic(0);
		A_xx += this->axis1(1) * this->axis1(1) * anisotropic(1);
		A_xx += this->axis1(2) * this->axis1(2) * anisotropic(2);

		// Get the lab-frame A_yy component
		double A_yy = this->axis2(0) * this->axis2(0) * anisotropic(0);
		A_yy += this->axis2(1) * this->axis2(1) * anisotropic(1);
		A_yy += this->axis2(2) * this->axis2(2) * anisotropic(2);

		return arma::cx_double((A_xx - A_yy) / 4.0, A_xy / 2.0);
	}
	// -----------------------------------------------------
	// Non-member non-friend functions
	// -----------------------------------------------------
	bool HasTrajectory(const Tensor &_tensor)
	{
		return (_tensor.TrajectoryLength() > 0);
	}

	bool IsIsotropic(const Tensor &_tensor)
	{
		// The tolerance used here is arbitrary. Use e.g. arma::eps or datum::eps for a better comparison
		if (std::abs(_tensor.Anisotropic()(0)) + std::abs(_tensor.Anisotropic()(1)) + std::abs(_tensor.Anisotropic()(2)) < 1e-10)
			return true;

		return false;
	}

	// -----------------------------------------------------
	// Non-member non-friend ActionTarget Check functions
	// -----------------------------------------------------
	// Make sure that the isotropic value has a valid value (not NaN or infinite)
	bool CheckActionScalarTensorIsotropicPart(const double &_d)
	{
		return std::isfinite(_d);
	}
	// -----------------------------------------------------
}
