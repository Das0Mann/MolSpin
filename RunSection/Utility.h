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
    
    typedef arma::cx_vec (*RungeKuttaFuncArma)(double t, arma::sp_cx_mat &, arma::cx_vec &, arma::cx_vec);
    
    /// Runge-Kutta-Fehlberg method (4th and 5th order) with adaptive time step control
    ///     @param L: Liouvillian superoperator (sparse complex matrix)
    ///     @param rho0: Initial density matrix (complex vector)
    ///     @param drhodt: Time derivative of the density matrix (complex vector)
    ///     @param dumpstep: Initial Time step (double)
    ///     @param func: Right hand side of the master equation (function pointer - see RungeKuttaFuncArma)
    ///     @param tolerance: Pair of tolerances for adaptive step size control (pair of doubles)
    ///     @param MinTimeStep: Minimum allowed time step (double) - Optional, default = 1e-6
    ///     @param MaxTimeStep: Maximum allowed time step (double) - Optional, default = 1e6
    ///     @param time: Current time (double) - Optional, default = 0
    ///     @return New time step (double)
    double RungeKutta45Armadillo(arma::sp_cx_mat &, arma::cx_vec &, arma::cx_vec &, double, RungeKuttaFuncArma, std::pair<double, double>, double MinTimeStep = 1e-6, double MaxTimeStep = 1e6, double time = 0);

    //DONT USE THESE FUNCTIONS THEY ARE SLOW 
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