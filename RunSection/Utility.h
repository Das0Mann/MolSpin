/////////////////////////////////////////////////////////////////////////
// Utility (RunSection module)
// ------------------
// Utility functions
// functions for multiple tasks
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////
#ifndef MOD_RunSection_Utility
#define MOD_RunSection_Utility

// #include <Eigen/Sparse>
// #include <Eigen/Core>
#include "SpinAPIfwd.h"

namespace RunSection
{
    // typedef Eigen::SparseMatrix<std::complex<double>, Eigen::RowMajor> Matrix;
    // typedef Eigen::Triplet<std::complex<double>, int32_t> T;

    // typedef Matrix(*RungeKuttaFunc)(Matrix&, Matrix&, Matrix);
    typedef arma::cx_vec (*RungeKuttaFuncArma)(double t, arma::sp_cx_mat &, arma::cx_vec &, arma::cx_vec);

    // Matrix ConvertArmadilloToEigen(arma::sp_cx_mat& ArmaMat);
    // Matrix ConvertAramdilloToEigen(arma::cx_vec& ArmaVec);
    //
    // arma::sp_cx_mat ConvertEigenToArmadillo(Matrix& EigenMat);
    //
    // double RungeKutta4AdaptiveTimeStepEigen(Matrix&, Matrix&, Matrix&, double, RungeKuttaFunc, std::pair<double, double>, double MinTimeStep = 1e-6);
    double RungeKutta45Armadillo(arma::sp_cx_mat &, arma::cx_vec &, arma::cx_vec &, double, RungeKuttaFuncArma, std::pair<double, double>, double MinTimeStep = 1e-6, double MaxTimeStep = 1e6, double time = 0);
}

#endif