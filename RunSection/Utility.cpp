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
                                                                          {12 / 13, {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0}},
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
}
