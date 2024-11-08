/////////////////////////////////////////////////////////////////////////
// Utility (RunSection module)
// ------------------
// Utility functions 
// functions for multiple tasks 
// 
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_Utility
#define MOD_RunSection_Utility

#include <Eigen/Sparse>
#include <Eigen/Core>
#include "SpinAPIfwd.h"

namespace RunSection
{
    typedef Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> Matrix;
    typedef Eigen::Triplet<std::complex<double>, int32_t> T;
    typedef Matrix(*RungeKuttaFunc)(Matrix&, Matrix&, Matrix);

    Matrix ConvertArmadilloToEigen(arma::sp_cx_mat& ArmaMat);
    arma::sp_cx_mat ConvertEigenToArmadillo(Matrix& EigenMat);

    double RungeKutta4AdaptiveTimeStepEigen(Matrix&, Matrix&, Matrix&, double, RungeKuttaFunc, std::pair<double, double>);
}

#endif