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

    //SparseMatrixSolvers
    //Preconditioned BiCGSTAB solver

    enum class PreconditionerType
    {
        None,
        IncompleteBiCGSTAB,
        SPAI,
        JACOBI,
        CUSTOM
    };

    arma::cx_vec BiCGSTAB(arma::sp_cx_mat &A, arma::cx_vec &b, PreconditionerType preconditioner, arma::sp_cx_mat K = arma::sp_cx_mat(), double tol = 1e-6, int max_iter = 1000, int max_preconditioner_iter = -1); //BiCGSTAB solver with preconditioner
    arma::sp_cx_mat IncompleteBiCGSTAB(arma::sp_cx_mat &A, int max_iter = 5); //Incomplete BiCGSTAB solver used to gererate a preconditioner
    arma::sp_cx_mat SPAI(arma::sp_cx_mat &A, int max_iter = 5); //Sparse Approximate Inverse preconditioner
    arma::sp_cx_mat JACOBI(arma::sp_cx_mat &A); //Jacobi preconditioner - leading  diagonal of A

    std::vector<int> LUDecomposition(arma::sp_cx_mat &K);
    arma::cx_vec LUSolve(arma::sp_cx_mat &K, std::vector<int> &P, arma::cx_vec &b); //LU decomposition and solve


}

#endif