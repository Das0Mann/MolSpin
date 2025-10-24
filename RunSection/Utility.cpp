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

    // Matrix ConvertArmadilloToEigen(arma::sp_cx_mat &ArmaMat)
    //{
    //     std::vector<T> Triplets(ArmaMat.n_nonzero);
    //
    //    int index = 0;
    //    #pragma omp parallel for schedule(static,1)
    //	for(unsigned int i = 0; i < ArmaMat.n_rows; i++)
    //	{
    //		arma::sp_cx_mat temp = ArmaMat.row(i);
    //		for(arma::sp_cx_mat::const_iterator it = temp.begin(); it != temp.end(); it++)
    //		{
    //            #pragma omp critical
    //            {
    //                Triplets[index] = T(i,it.col(),(*it));
    //                index++;
    //            }
    //		}
    //	}
    //
    //    Matrix ReturnMat(ArmaMat.n_rows, ArmaMat.n_rows);
    //	ReturnMat.setFromTriplets(Triplets.begin(), Triplets.end());
    //    return ReturnMat;
    //}

    // Matrix ConvertAramdilloToEigen(arma::cx_vec &ArmaVec)
    //{
    //     std::vector<T> Triplets(ArmaVec.n_elem);
    //     int index = 0;
    //
    //    #pragma omp parallel for schedule(static,1)
    //    for(unsigned int i = 0; i < ArmaVec.n_rows; i++)
    //    {
    //        Triplets[i] = T(i,0,ArmaVec[i]);
    //    }
    //
    //    Matrix ReturnMat(ArmaVec.n_rows, 1);
    //    ReturnMat.setFromTriplets(Triplets.begin(), Triplets.end());
    //    return ReturnMat;
    //}
    //
    // arma::sp_cx_mat ConvertEigenToArmadillo(Matrix &EigenMat)
    //{
    //    unsigned int a = EigenMat.nonZeros();
    //    arma::umat positions(2,a);
    //    arma::cx_vec values(a);
    //    int ind = 0;
    //    #pragma omp parallel for schedule(static,1)
    //    for(unsigned int i = 0; i < EigenMat.rows(); i++)
    //    {
    //        Eigen::SparseVector<std::complex<double>> vec = EigenMat.row(i);
    //        unsigned int s = vec.nonZeros();
    //        for(unsigned int e = 0; e < s; e++)
    //        {
    //            int index = vec.data().index(e);
    //            std::complex val = vec.coeff(index);
    //            #pragma omp critical
    //            {
    //                positions(0, ind) = i;
    //                positions(1,ind) = e;
    //                values(ind) = val;
    //                ind++;
    //            }
    //        }
    //    }
    //
    //    return arma::sp_cx_mat(positions, values);
    //}
    //
    // double RungeKutta4AdaptiveTimeStepEigen(Matrix& L, Matrix& rho0, Matrix& drhodt, double dumpstep, RungeKuttaFunc func, std::pair<double, double> tolerance, double MinTimeStep)
    //{
    //    std::vector<double> timesteps = {dumpstep*0.5, dumpstep, dumpstep * 2};
    //    Matrix k0(rho0.rows(),1);
    //
    //    auto RungeKutta4 = [&k0](Matrix& L1, Matrix& rho01, double timestep1, RungeKuttaFunc func1) {
    //        Matrix k1(rho01.rows(),1);
    //        Matrix k2(rho01.rows(),1);
    //        Matrix k3(rho01.rows(),1);
    //        Matrix k4(rho01.rows(),1);
    //
    //        k1 = func1(L1, k0, rho01);
    //        {
    //            Matrix temp = std::complex<double>(0.5*timestep1,0) * k1;
    //            k2 = func1(L1, temp, rho01);
    //        }
    //        {
    //            Matrix temp = std::complex<double>(0.5*timestep1,0) * k2;
    //            k3 = func1(L1, temp, rho01);
    //        }
    //        {
    //            Matrix temp = std::complex<double>(timestep1,0) * k3;
    //            k4 = func1(L1, temp, rho01);
    //        }
    //
    //        {
    //            Matrix temp1 = std::complex<double>(2.0,0.0) * k2;
    //            Matrix temp2 = std::complex<double>(2.0,0.0) * k3;
    //            Matrix temp3 = k1 + temp1 + temp2 + k4;
    //            Matrix temp4 = (timestep1/6.0) * temp3;
    //            Matrix temp5 = rho01+temp4;
    //            return temp5;
    //        }
    //        //return rho0 + ((timestep/6.0) * (k1 + (std::complex<double>(2.0,0.0) * k2) + (std::complex<double>(2.0,0.0) * k3) + k4));
    //    };
    //
    //    auto Change = [&rho0](const Matrix rho1) {
    //        Matrix diff = rho0 - rho1;
    //        double abs = std::abs(diff.sum());
    //        double abs2 = std::abs(rho0.sum());
    //        double change = abs/abs2;
    //        return change;
    //    };
    //
    //    std::vector<Matrix> steps(3);
    //    if(timesteps[0] > MinTimeStep)
    //    {
    //        #pragma omp parallel for schedule(static,1)
    //        for (int i = 0; i < 3; i++)
    //        {
    //            steps[i] = RungeKutta4(L, rho0, timesteps[i], func);
    //        }
    //    }
    //    else
    //    {
    //        #pragma omp parallel for schedule(static,1)
    //        for (int i = 1; i < 3; i++)
    //        {
    //            steps[i] = RungeKutta4(L, rho0, timesteps[i], func);
    //        }
    //    }
    //
    //    if(timesteps[0] > MinTimeStep)
    //    {
    //        if(Change(steps[0]) > tolerance.second)
    //        {
    //            drhodt = steps[0];
    //            return timesteps[0];
    //        }
    //    }
    //    else if(Change(steps[2]) < tolerance.first)
    //    {
    //        drhodt = steps[2];
    //        return timesteps[2];
    //    }
    //    else
    //    {
    //        drhodt = steps[1];
    //        return dumpstep;
    //    }
    //
    //    drhodt = RungeKutta4(L, rho0, timesteps[1], func);
    //
    //    return dumpstep;
    //}

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

    arma::cx_vec ThomasBlockSolver(arma::sp_cx_mat &A, arma::cx_vec &b, int block_size)
    {
        arma::cx_vec x = arma::cx_vec(arma::size(b), arma::fill::zeros);
        return x;
    }

    bool BlockSolver(arma::sp_cx_mat &A, arma::cx_vec &b, int block_size, arma::cx_vec &x)
    {
        bool inverted = false;
        if(IsBlockTridiagonal(A))
        {
            //return ThomasBlockSolver(A, b, block_size);
            std::cout << "Block Thomas Solver not yet implemented, using Block Matrix Inverse instead." << std::endl;
        }

        arma::cx_mat A_inv = BlockMatrixInverse(A, block_size, inverted);
        if(!inverted)
        {
            x = arma::cx_vec(arma::size(b), arma::fill::zeros);
            return false;
        }
        x = A_inv * b;
        return true;
    }

    arma::cx_mat BlockMatrixInverse(arma::sp_cx_mat &A, int block_size, bool &Invertible)
    {
        //Matrix Partitions
        arma::sp_cx_mat A11, A12, A21, A22;
        A11 = A.submat(0, 0, block_size -1, block_size -1);
        A12 = A.submat(0, block_size, block_size -1, A.n_cols -1);
        A21 = A.submat(block_size, 0, A.n_rows-1, block_size -1);
        A22 = A.submat(block_size, block_size, A.n_rows -1, A.n_cols -1);

        if (A22.n_rows > block_size)
        {
            BlockMatrixInverse(A22, block_size, Invertible);
        }

        //Check if A11 and A22 are invertible
        arma::cx_mat A11_inv, A22_inv; //The inverse of a sparse matrix is usually dense, so we use a dense matrix here
        //Check invertibility wihtin a scope, that way if not invertable we don't keep the failed inverse
        bool A11_invertible, A22_invertible;
        bool SchurComplementA, SchurComplementB, BothComplements;
        {
            bool A11_invertible = arma::inv(A11_inv, arma::cx_mat(A11));
            bool A22_invertible = arma::inv(A22_inv, arma::cx_mat(A22));

            if(!A11_invertible)
            {
                A11_inv = arma::cx_mat(); 
            }
            if(!A22_invertible)
            {
                A22_inv = arma::cx_mat(); 
            }
        }

        bool SchurComplementA = A11_invertible && !A22_invertible;
        bool SchurComplementB = !A11_invertible && A22_invertible;
        bool BothComplements = A11_invertible && A22_invertible;
        Invertible = true;

        if(SchurComplementA)
        {
            return SchurComplementA(A11_inv, A12, A21, A22, Invertible);
        }
        else if(SchurComplementB)
        {
            return SchurComplementB(A11, A12, A21, A22_inv, Invertible);
        }
        else if(BothComplements)
        {
            return BothSchurComponents(A11_inv, A12, A21, A22_inv, Invertible);
        }
        else
        {
            Invertible = false;
            return arma::cx_mat();
        }

    }

    arma::cx_mat SchurComplementA(arma::cx_mat &A11_inv, arma::sp_cx_mat &A12, arma::sp_cx_mat &A21, arma::sp_cx_mat &A22, bool &Invertible)
    {
        arma::cx_mat S = A22 - A21 * A11_inv * A12;
        arma::cx_mat S_inv;
        bool S_invertible = arma::inv(S_inv, S);
        if(!S_invertible)
        {
            Invertible = false;
            return arma::cx_mat();
        }
        //Construct the inverse matrix using the Schur complement
        arma::cx_mat P11 = A11_inv + A11_inv * A12 * S_inv * A21 * A11_inv;
        arma::cx_mat P12 = -1 * A11_inv * A12 * S_inv;
        arma::cx_mat P21 = -1 * S_inv * A21 * A11_inv;
        arma::cx_mat P22 = S_inv;

        arma::cx_mat Inv = arma::cx_mat(A11_inv.n_rows * 2, A11_inv.n_cols * 2);
        Inv.submat(0, 0, A11_inv.n_rows -1, A11_inv.n_cols -1) = P11;
        Inv.submat(0, A11_inv.n_cols, A11_inv.n_rows -1, Inv.n_cols -1) = P12;
        Inv.submat(A11_inv.n_rows, 0, Inv.n_rows -1, A11_inv.n_cols -1) = P21;
        Inv.submat(A11_inv.n_rows, A11_inv.n_cols, Inv.n_rows -1, Inv.n_cols -1) = P22;

        Invertible = true;
        return Inv;
    }

    arma::cx_mat SchurComplementB(arma::sp_cx_mat &A11, arma::sp_cx_mat &A12, arma::cx_mat &A21, arma::cx_mat &A22_inv, bool &Invertible)
    {
        arma::cx_mat S = A11 - A12 * A22_inv * A21;
        arma::cx_mat S_inv;
        bool S_invertible = arma::inv(S_inv, S);
        if(!S_invertible)
        {
            Invertible = false;
            return arma::cx_mat();
        }
        //Construct the inverse matrix using the Schur complement
        arma::cx_mat P11 = S_inv;
        arma::cx_mat P12 = -1 * S_inv * A12 * A22_inv;
        arma::cx_mat P21 = -1 * A22_inv * A21 * S_inv;
        arma::cx_mat P22 = A22_inv + A22_inv * A21 * S_inv * A12 * A22_inv;

        arma::cx_mat Inv = arma::cx_mat(A11.n_rows * 2, A11.n_cols * 2);
        Inv.submat(0, 0, A11.n_rows -1, A11.n_cols -1) = P11;
        Inv.submat(0, A11.n_cols, A11.n_rows -1, Inv.n_cols -1) = P12;
        Inv.submat(A11.n_rows, 0, Inv.n_rows -1, A11.n_cols -1) = P21;
        Inv.submat(A11.n_rows, A11.n_cols, Inv.n_rows -1, Inv.n_cols -1) = P22;

        bool Invertible = true;
        return Inv;
    }

    arma::cx_mat BothSchurComponents(arma::cx_mat&A11, arma::cx_mat &A11_inv, arma::sp_cx_mat &A12, arma::sp_cx_mat &A21, arma::cx_mat &A22, arma::cx_mat &A22_inv, bool &Invertible)
    {
        arma::cx_mat S1 = A11 - A12 * A22_inv * A21;
        arma::cx_mat S2 = A22 - A21 * A11_inv * A12;
        arma::cx_mat S1_inv, S2_inv;
        bool S1_invertible = arma::inv(S1_inv, S1);
        bool S2_invertible = arma::inv(S2_inv, S2);
        if(!S1_invertible || !S2_invertible)
        {
            Invertible = false;
            return arma::cx_mat();
        }
        //Construct the inverse matrix using the Schur complement
        arma::cx_mat P11 = S1_inv;
        arma::cx_mat P12 = -1 * S1_inv * A12 * A22_inv;
        arma::cx_mat P21 = -1 * S2_inv * A21 * A11_inv;
        arma::cx_mat P22 = S2_inv;

        arma::cx_mat Inv = arma::cx_mat(A11.n_rows * 2, A11.n_cols * 2);
        Inv.submat(0, 0, A11.n_rows -1, A11.n_cols -1) = P11;
        Inv.submat(0, A11.n_cols, A11.n_rows -1, Inv.n_cols -1) = P12;
        Inv.submat(A11.n_rows, 0, Inv.n_rows -1, A11.n_cols -1) = P21;
        Inv.submat(A11.n_rows, A11.n_cols, Inv.n_rows -1, Inv.n_cols -1) = P22;
        
        bool Invertible = true;
        return Inv;
    }

    bool IsBlockTridiagonal(arma::sp_cx_mat &A, int block_size)
    {
        int n_blocks = A.n_rows / block_size;
        if n_blocks == 2
        {
            return true;   
        }

        for (int row = 0; row < n_blocks; row++)
        {
            for (int col = 0; col < n_blocks; col++)
            {
                if(std::abs(row-col) > 1)
                {
                    arma::sp_cx_mat block = A.submat(row*block_size, col*block_size, (row+1)*block_size -1, (col+1)*block_size -1);
                    arma::vec non_zero = arma::nonzeros(block); //should be very efficient for sparse matrices
                    if(non_zero.n_elem > 0)
                    {
                        return false;
                    }
                }
            }
        }

        return true;
    }

    arma::cx_vec BiCGSTAB(arma::sp_cx_mat &A, arma::cx_vec &b, PreconditionerType preconditoner ,arma::sp_cx_mat K, double tol, int max_iter, int max_preconditoner_iter)
    {

        if(preconditoner == PreconditionerType::None)
        {
            K = arma::sp_cx_mat(arma::size(A));
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


        arma::sp_cx_mat K_1, K_2;
        auto P = LUDecomposition(K);
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
            arma::cx_vec y = LUSolve(K,P, rho_prev);
            arma::cx_vec v = A * y;
            arma::cx_double alpha = rho_k_1 / arma::dot(r_naught_hat, v);
            arma::cx_vec h = x + alpha * y;
            arma::cx_vec s = r - alpha * v;
            if(arma::norm(s) < tol)
            {
                std::cout << "Converged in " << k << " iterations." << std::endl;
                return x;
            }
            arma::cx_vec z = LUSolve(K,P,s);
            arma::cx_vec t = A * z;
            arma::cx_double omega = arma::dot(t, s) / arma::dot(t, t);
            x = h + omega * z;
            r = s - omega * t;
            if(arma::norm(r) < tol)
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
            x.col(i) = b;
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
