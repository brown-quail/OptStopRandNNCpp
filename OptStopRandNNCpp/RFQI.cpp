#include "RFQI.h"
#include <Eigen/Dense>
#include <math.h>
#include <iostream>


using namespace std;
// do RFQI. first input: payoff function pointer
double doRFQI(payoff g, double strike, 
	double*** path, int m, int d, int N, 
	double discount, int K, double T,
	int max_iter, double tol, mt19937_64 urng)
{
	double p0 = 0;

	// randomized NN fix
	Eigen::MatrixXd A = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, d + 2, urng);
	Eigen::MatrixXd b = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, 1, urng);

	// preparation
	Eigen::MatrixXd theta(K, 1);
	Eigen::MatrixXd next_theta(K, 1);
	theta.setZero();
	next_theta.setZero();

	double** pni = new double* [N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		pni[i] = new double[2 * m];
	}

	Eigen::MatrixXd phi(K, 1);
	Eigen::MatrixXd phi_next(K, 1);
	Eigen::MatrixXd xtild(d + 2, 1);
	Eigen::MatrixXd xtild_next(d + 2, 1);
	Eigen::MatrixXd x(d, 1);
	Eigen::MatrixXd phiphiT(K, K);
	Eigen::MatrixXd phipni(K, 1);
	int iter = 0;

	double alpha = exp(-discount * T / N);

	bool inverseFlag = false;
	Eigen::MatrixXd phiphiTinverse;
	do
	{
		theta = next_theta;
		phiphiT.setConstant(0);
		phipni.setConstant(0);
		for (int i = 0; i < 2 * m; i++)
		{
			pni[N][i] = g(path, i, N, d, strike);
			for (int n = N - 1; n > 0; n--)
			{
				for (int j = 0; j < d; j++)
				{
					xtild(j, 0) = path[i][j][n];
				}
				xtild(d, 0) = (double)n;
				xtild(d + 1, 0) = ((double)N - (double)n);
				Eigen::MatrixXd temp(K - 1, 1);
				temp = leakyRelu(A * xtild + b);
				temp.conservativeResize(K, 1);
				temp(K - 1, 0) = 1;
				phi = temp;
				
				pni[n][i] = max(g(path, i, n, d, strike), (phi.transpose() * theta)(0, 0));
				//pni[n][i] = g(path, i, n, d, strike) >= (theta.transpose() * phi)(0, 0) ? g(path, i, n, d, strike) : alpha * pni[n + 1][i];
				
				if (i < m)
				{
					phiphiT = phiphiT + phi * phi.transpose();
					phipni = phipni + phi * pni[n + 1][i];
				}
			}
		}

		if (inverseFlag == false)
		{
			phiphiTinverse = phiphiT.inverse();
			inverseFlag = true;
		}
		next_theta = alpha * phiphiTinverse * phipni;
		
		//next_theta = phiphiT.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(alpha * phipni);
#ifdef _DEBUG
		std::cout << "[iteration " << iter << "] ";
		std::cout << "theta difference: " << (theta - next_theta).norm() << endl;
#endif
	} while (iter++ < max_iter && (theta - next_theta).norm() > tol);

	/*
	for (int i = 0; i < m; i++)
	{
		p0 += pni[1][i + m];
	}
	p0 /= m;

	p0 = max(g(path, 0, 0, d, strike), p0);
	*/
		
	double* price = new double[m];
	for (int i = 0; i < m; i++)
	{
		price[i] = 0;
	}
	for (int i = 0; i < m; i++)
	{
		for (int n = 0; n < N + 1; n++)
		{
			for (int j = 0; j < d; j++)
			{
				xtild(j, 0) = path[i + m][j][n];
			}
			xtild(d, 0) = (double)n;
			xtild(d + 1, 0) = ((double)N - (double)n);
			Eigen::MatrixXd temp(K - 1, 1);
			temp = leakyRelu(A * xtild + b);
			temp.conservativeResize(K, 1);
			temp(K - 1, 0) = 1;
			phi = temp;

			double pay = g(path, i + m, n, d, strike);
			double conti_value = (theta.transpose() * phi)(0, 0);
			conti_value = std::max(conti_value, 0.0);
			if (pay > conti_value || (n == N))
			{
				price[i] = pay;
				for (int k = 0; k < n; k++)
				{
					price[i] = price[i] * alpha;
				}
				break;
			}
		}
	}
	p0 = 0;
	for (int i = 0; i < m; i++)
	{
		p0 += price[i];
	}
	p0 /= m;
	


	for (int k = 0; k < N + 1; k++)
	{
		delete[] pni[k];
	}
	delete[] pni;
	delete[] price;

	return p0;
}