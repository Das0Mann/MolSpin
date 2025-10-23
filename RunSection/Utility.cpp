/////////////////////////////////////////////////////////////////////////
// Utility implementation (RunSection module)
//
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2025 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////

#include "Utility.h"
namespace RunSection
{
    typedef arma::sp_cx_mat MatrixArma;
    typedef arma::cx_vec VecType;

    double RungeKutta45Armadillo(arma::sp_cx_mat &L, arma::cx_vec &rho0, arma::cx_vec &drhodt, double dumpstep, RungeKuttaFuncArma func, std::pair<double, double> tolerance, double MinTimeStep, double MaxTimeStep, double time)
    {
        VecType k0(rho0.n_rows);

        std::vector<std::pair<float, std::vector<float>>> ButcherTable = {{0.0, {}},
                                                                          {0.25, {0.25}},
                                                                          {3.0 / 8.0, {3.0 / 32.0, 9.0 / 32.0}},
                                                                          {12.0 / 13.0, {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0}},
                                                                          {1.0, {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0}},
                                                                          {1.0 / 2.0, {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0}},
                                                                          {0.0, {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0}},
                                                                          {0.0, {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0, 0.0}}};
        auto RungeKutta45 = [&k0, &ButcherTable, &time](MatrixArma &L1, VecType &rho01, double t, RungeKuttaFuncArma func1)
        {
            VecType k1(rho01.n_rows);
            VecType k2(rho01.n_rows);
            VecType k3(rho01.n_rows);
            VecType k4(rho01.n_rows);
            VecType k5(rho01.n_rows);
            VecType k6(rho01.n_rows);

            std::vector<VecType> kvec = {k1, k2, k3, k4, k5, k6};

            auto GetK = [&ButcherTable](int index, std::vector<VecType> kv)
            {
                VecType temp(kv[0].n_rows);
                for (int e = 0; e < int(ButcherTable[index].second.size()); e++)
                {
                    temp = temp + (ButcherTable[index].second[e] * kv[e]);
                }
                return temp;
            };

            int i = 0;
            kvec[0] = t * func1(time + ButcherTable[i].first, L1, k0, rho01);
            // std::cout << kvec[0] << std::endl;
            i += 1;
            for (; i < 6; i++)
            {
                VecType temp = GetK(i, kvec);
                // std::cout << temp << std::endl;
                kvec[i] = t * func1(time + ButcherTable[i].first, L1, temp, rho01);
            }

            VecType ReturnVecRK4 = rho01;
            for (i = 0; i < int(ButcherTable[7].second.size()); i++)
            {
                // std::cout << i << std::endl;
                ReturnVecRK4 += (ButcherTable[7].second[i] * kvec[i]);
            }

            VecType ReturnVecRK5 = rho01;
            for (i = 0; i < int(ButcherTable[6].second.size()); i++)
            {
                ReturnVecRK5 += (ButcherTable[6].second[i] * kvec[i]);
            }

            return std::make_tuple(ReturnVecRK4, ReturnVecRK5);
        };

        auto [RK4, RK5] = RungeKutta45(L, rho0, dumpstep, func);

        double change = 0;
        {
            VecType diff = RK5 - RK4;
            double sum = 0;

#pragma omp parallel for reduction(+ : sum)
            for (int i = 0; i < int(diff.n_rows); i++)
            {
                sum += std::pow(std::abs(diff[i]), 2);
            }

            change = std::sqrt(sum);
        }

        auto Adjusth = [](double tol, double ch)
        {
            double h4 = (tol / (2 * ch));
            return std::sqrt(std::sqrt(h4));
        };

        double NewStepSize = 0.0;
        if (change < tolerance.first && dumpstep < MaxTimeStep)
        {
            NewStepSize = dumpstep * Adjusth(tolerance.first, change);
            if (NewStepSize > MaxTimeStep)
            {
                NewStepSize = MaxTimeStep;
            }
        }
        else if (change > tolerance.second && dumpstep > MinTimeStep)
        {
            NewStepSize = dumpstep * Adjusth(tolerance.second, change);
            if (NewStepSize < MinTimeStep)
            {
                NewStepSize = MinTimeStep;
            }
        }
        else
        {
            NewStepSize = dumpstep;
        }

        drhodt = RK4;
        return NewStepSize;
    }


    //DONT USE THESE FUNCTIONS THEY ARE SLOW
    arma::cx_vec BiCGSTAB(arma::sp_cx_mat &A, arma::cx_vec &b, PreconditionerType preconditoner ,arma::sp_cx_mat K, double tol, int max_iter, int max_preconditoner_iter)
    {

        if(preconditoner == PreconditionerType::None)
        {
            K = arma::sp_cx_mat(arma::size(A));
            K.eye(); // No preconditioner, use identity matrix
        }
        else if(preconditoner == PreconditionerType::IncompleteBiCGSTAB)
        {
            if (max_preconditoner_iter < 0)
            {
                max_preconditoner_iter = 5; 
            }
            K = IncompleteBiCGSTAB(A, max_preconditoner_iter);
        }
        else if(preconditoner == PreconditionerType::SPAI)
        {
            if (max_preconditoner_iter < 0)
            {
                max_preconditoner_iter = 50; 
            }
            
            K = SPAI(A, max_preconditoner_iter);
        }
        else if(preconditoner == PreconditionerType::JACOBI)
        {
            K = JACOBI(A);
        }


        arma::cx_mat K_1, K_2;
        //auto P = LUDecomposition(K);
        arma::cx_mat DenseK = arma::conv_to<arma::cx_mat>::from(K);
        arma::lu(K_1, K_2, DenseK);
        //arma::lu(K_1, K_2, K);
        //arma::lu()

        arma::cx_vec x = arma::cx_vec(arma::size(b), arma::fill::zeros);
        arma::cx_vec r_naught = b - A * x;
        arma::cx_vec r_naught_hat = r_naught;
        arma::cx_vec r = r_naught;
        arma::cx_vec rho_naught = r_naught;
        arma::cx_vec rho_prev = rho_naught;
        arma::cx_double rho_k_1 = arma::dot(rho_naught, rho_naught);

        for(int k = 1; k <= max_iter; k++)
        {
            //arma::cx_vec y = LUSolve(K,P, rho_prev); //too slow
            arma::cx_vec y = arma::cx_vec(arma::size(rho_prev), arma::fill::zeros);
            arma::solve(y,DenseK, rho_prev);
            arma::cx_vec v = A * y;
            arma::cx_double alpha = rho_k_1 / arma::dot(r_naught_hat, v);
            arma::cx_vec h = x + alpha * y;
            arma::cx_vec s = r - alpha * v;
            double norms = arma::norm(s);
            if(norms < tol)
            {
                std::cout << "Converged in " << k << " iterations." << std::endl;
                return x;
            }
            //arma::cx_vec z = LUSolve(K,P,s);
            arma::cx_vec z = arma::cx_vec(arma::size(s), arma::fill::zeros);
            arma::solve(z, DenseK, s);
            arma::cx_vec t = A * z;
            arma::cx_vec K_1_invt = arma::cx_vec(arma::size(t), arma::fill::zeros);
            arma::solve(K_1_invt, K_1, t);
            arma::cx_vec K_1_invs = arma::cx_vec(arma::size(s), arma::fill::zeros);
            arma::solve(K_1_invs, K_1, s);
            arma::cx_double omega = arma::dot(K_1_invt, K_1_invs) / arma::dot(K_1_invt, K_1_invt);
            x = h + omega * z;
            r = s - omega * t;
            double normr = arma::norm(r);
            if(normr < tol)
            {
                std::cout << "Converged in " << k << " iterations." << std::endl;
                return x;
            }
            arma::cx_double rho_k = arma::dot(r_naught_hat, r);
            arma::cx_double beta = (rho_k / rho_k_1) * (alpha / omega);
            rho_prev = r + beta * (rho_prev - omega * v);
            rho_k_1 = rho_k;
        }
        std::cout << "Did not converge in " << max_iter << " iterations." << std::endl;
        return x; // Return the last computed x, even if it did not converge
    }

    arma::sp_cx_mat IncompleteBiCGSTAB(arma::sp_cx_mat &A, int max_iter)
    {
        int n = A.n_rows;
        arma::sp_cx_mat I = arma::sp_cx_mat(arma::size(A));
        I.eye();
        arma::sp_cx_mat x = arma::sp_cx_mat(arma::size(A));
        for(int i = 0; i < n; i++)
        {
            arma::cx_vec col = arma::cx_vec(n, arma::fill::zeros);
            for (int j = 0; j < n; j++)
            {
                col(j) = A(i,j);
            }
            arma::cx_vec b = BiCGSTAB(A, col, PreconditionerType::None, arma::sp_cx_mat(), max_iter = max_iter);
            {
                x.col(i) = b;
            }
        }
        return x;
    }

    arma::sp_cx_mat SPAI(arma::sp_cx_mat &A, int max_iter)
    {
        arma::sp_cx_mat I = arma::sp_cx_mat(arma::size(A));
        I.eye();
        arma::cx_double alpha = 2.0 / arma::norm(A * arma::trans(A), 1);
        arma::sp_cx_mat M = alpha * A;
        for(int i = 0; i < max_iter; i++)
        {
            arma::sp_cx_mat C = A * M;
            arma::sp_cx_mat G = I - C;
            arma::sp_cx_mat AG = A * G;
            arma::cx_double trace = arma::trace(arma::trans(G) * AG);
            arma::cx_double norm = arma::norm(AG, 1);
            alpha = trace / std::pow(norm,2);
            M = M + alpha * G;
        }
        return M;
    }

    arma::sp_cx_mat JACOBI(arma::sp_cx_mat &A)
    {
        arma::sp_cx_mat K = arma::sp_cx_mat(arma::size(A));
        K = A.diag();
        return K;
    }

    std::vector<int> LUDecomposition(arma::sp_cx_mat &A)
    {
        arma::sp_cx_mat L, U;
        int n = A.n_rows;
        std::vector<int> permuation;
        for(int i = 0; i <= n; i++)
        {
            permuation.push_back(i);
        }

        for (int i = 0; i < n; i++)
        {
            double max_val = 0.0;
            int max_index = i;

            for (int k = i; k < n; k++)
            {
                std::complex<double> ki = A(k,i);
                double val = std::abs(ki);
                if (val > max_val)
                {
                    max_val = val;
                    max_index = k;
                }
            }

            if (max_index != i)
            {
                int j = permuation[i];
                permuation[i] = permuation[max_index];
                permuation[max_index] = j;
                A.swap_rows(i, max_index);
                permuation[n] = permuation[n] + 1; // Increment the permutation count
            }

            for (int j = i + 1; j < n; j++)
            {
                std::complex<double> ii,ji,ik,jk;
                ii = A(i,i);
                ji = A(j,i);
                ik = A(i,j);
                if (ii == std::complex<double>(0, 0))
                {
                    throw std::runtime_error("Matrix is singular, cannot perform LU decomposition.");
                }
                A(j,i) = ji / ii;
                for (int k = i + 1; k < n; k++)
                {
                    jk = A(j,k);
                    A(j,k) = jk - (ji * ik);
                }
            }
        }
        return permuation;
    }

    arma::cx_vec LUSolve(arma::sp_cx_mat &K, std::vector<int> &P, arma::cx_vec &b)
    {
        int n = K.n_rows;
        arma::cx_vec x = arma::cx_vec(n, arma::fill::zeros);
        for (int i = 0; i < n; i++)
        {
            x(i) = b(P[i]);
            for (int j = 0; j < i; j++)
            {
                std::complex<double> xi, xj,ij;
                xi = x(i);
                xj = x(j);
                ij = K(i,j);
                x(i) = xi - ij * xj;
            }
        }

        for (int i = n - 1; i >= 0; i--)
        {
            for(int j = i + 1; j < n; j++)
            {
                std::complex<double> xi, xj, ij;
                xi = x(i);
                xj = x(j);
                ij = K(i,j);
                x(i) = x(i) - ij * xj;
            }
            std::complex<double> xi, ii;
            xi = x(i);
            ii = K(i,i);
            x(i) = x(i) / ii;
        }
        return x;
    }
}
