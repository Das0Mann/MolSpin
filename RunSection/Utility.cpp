/////////////////////////////////////////////////////////////////////////
// Utility implementation (RunSection module)
//
//
// Molecular Spin Dynamics Software - developed by Claus Nielsen and Luca Gerhards.
// (c) 2019 Quantum Biology and Computational Physics Group.
// See LICENSE.txt for license information.
/////////////////////////////////////////////////////////////////////////

#include "Utility.h"
namespace RunSection
{
    Matrix ConvertArmadilloToEigen(arma::sp_cx_mat &ArmaMat)
    {
        std::vector<T> Triplets(ArmaMat.n_nonzero);

        int index = 0;
        #pragma omp parallel for schedule(static,1)
    	for(unsigned int i = 0; i < ArmaMat.n_rows; i++)
    	{
    		arma::sp_cx_mat temp = ArmaMat.row(i);
    		for(arma::sp_cx_mat::const_iterator it = temp.begin(); it != temp.end(); it++)
    		{
                #pragma omp critical
                {
                    Triplets[index] = T(i,it.col(),(*it));
                    index++;
                }
    		}
    	}

        Matrix ReturnMat(ArmaMat.n_rows, ArmaMat.n_rows);
    	ReturnMat.setFromTriplets(Triplets.begin(), Triplets.end());
        return ReturnMat;
    }

    arma::sp_cx_mat ConvertEigenToArmadillo(Matrix &EigenMat)
    {
        unsigned int a = EigenMat.nonZeros();
        arma::umat positions(2,a);
        arma::cx_vec values(a);
        int ind = 0;
        #pragma omp parallel for schedule(static,1)
        for(unsigned int i = 0; i < EigenMat.rows(); i++)
        {
            Eigen::SparseVector<std::complex<double>> vec = EigenMat.row(i);
            unsigned int s = vec.nonZeros();
            for(unsigned int e = 0; e < s; e++)
            {
                int index = vec.data().index(e);
                std::complex val = vec.coeff(index);
                #pragma omp critical
                {
                    positions(0, ind) = i;
                    positions(1,ind) = e;
                    values(ind) = val;
                    ind++;
                }
            }
        }

        return arma::sp_cx_mat(positions, values);
    }

    double RungeKutta4AdaptiveTimeStepEigen(Matrix& L, Matrix& rho0, Matrix& drhodt, double dumpstep, RungeKuttaFunc func, std::pair<double, double> tolerance)
    {
        std::vector<double> timesteps = {dumpstep*0.5, dumpstep, dumpstep * 2};
        Matrix k0(rho0.rows(),1);

        auto RungeKutta4 = [&k0](Matrix& L, Matrix& rho0, double timestep, RungeKuttaFunc func) {
            Matrix k1(rho0.rows(),1);
            Matrix k2(rho0.rows(),1);
            Matrix k3(rho0.rows(),1);
            Matrix k4(rho0.rows(),1);

            k1 = func(L, k0, rho0);
            {
                Matrix temp = std::complex<double>(0.5*timestep,0) * k1;
                k2 = func(L, temp, rho0);
            }  
            {
                Matrix temp = std::complex<double>(0.5*timestep,0) * k2;
                k3 = func(L, temp, rho0);
            }
            {
                Matrix temp = std::complex<double>(timestep,0) * k3;
                k4 = func(L, temp, rho0);
            }

            return rho0 + ((timestep/6.0) * (k1 + (std::complex<double>(2.0,0.0) * k2) + (std::complex<double>(2.0,0.0) * k3) + k4));
        };

        //auto Change = [&rho0](const Matrix rho1) {
        //    Matrix diff = rho 0- rho1;
        //    double abs = diff.cwiseAbs();
        //    return (abs/rho0.cwiseAbs());
        //};

        //std::vector<Matrix> steps(3); 
        //for (int i = 0; i < 3; i++)
        //{
        //    steps[i] = RungeKutta4(L, rho0, timesteps[i], func);
        //}
//
        //if(Change(steps[0]) > tolerance.second)
        //{
        //    drhodt = steps[0];
        //    return timesteps[0];
        //}
        //else if(Change(steps[0]) < tolerance.first)
        //{
        //    drhodt = steps[2];
        //    return timesteps[2];
        //}
        //else
        //{
        //    drhodt = steps[1];
        //    return dumpstep;
        //}

        return dumpstep;
    }
}