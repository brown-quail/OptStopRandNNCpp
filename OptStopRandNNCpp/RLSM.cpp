#include "RLSM.h"
#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <Eigen/SVD>


double doRLSM(payoff g, double strike, 
	double*** path, int m, int d, int N,
	double discount, int K, double T,
	mt19937_64 urng)
{
	Eigen::MatrixXd A = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, d, urng);
	Eigen::MatrixXd b = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, 1, urng);

	Eigen::MatrixXd theta(K, 1);

	theta.setZero();

	double** pni = new double* [N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		pni[i] = new double[2 * m];
	}

	Eigen::MatrixXd phi(K, 1);
	Eigen::MatrixXd phiphiT(K, K);
	Eigen::MatrixXd phipni(K, 1);
	Eigen::MatrixXd x(d, 1);

	Eigen::MatrixXd phiM(K, m);
	Eigen::MatrixXd pniM(m, 1);

	double alpha = exp(-discount * T / N);

	for (int i = 0; i < 2 * m; i++)
	{
		pni[N][i] = g(path, i, N, d, strike);
	}
	for (int n = N - 1; n > 0; n--)
	{
		phiphiT.setConstant(0);
		phipni.setConstant(0);
		phiM.setConstant(0);
		for (int i = 0; i < 2 * m; i++)
		{
			for (int j = 0; j < d; j++)
			{
				x(j, 0) = path[i][j][n];
			}
			Eigen::MatrixXd temp(K - 1, 1);
			temp = leakyRelu(A * x + b);
			temp.conservativeResize(K, 1);
			temp(K - 1, 0) = 1;
			phi = temp;

			if (i < m)
			{
				phiM.col(i) = temp;
				pniM(i) = pni[n + 1][i];
			}
		}
		Eigen::MatrixXd phiMT = phiM.transpose();	// m x K
		// solve phiMT * theta = alpha * pniM
		// dimension (m x K) * (K, 1) = (m, 1)
		//theta = phiMT.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(alpha * pniM);
		theta = phiMT.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(alpha * pniM);
		for (int i = 0; i < 2 * m; i++)
		{
			pni[n][i] = g(path, i, n, d, strike) >= (theta.transpose() * phi)(0, 0) ? g(path, i, n, d, strike) : alpha * pni[n + 1][i];
		}
	}

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
				x(j, 0) = path[i + m][j][n];
			}
			Eigen::MatrixXd temp(K - 1, 1);
			temp = leakyRelu(A * x + b);
			temp.conservativeResize(K, 1);
			temp(K - 1, 0) = 1;
			phi = temp;

			double payoff = g(path, i + m, n, d, strike);
			double conti_value = (theta.transpose() * phi)(0, 0);
			conti_value = std::max(conti_value, 0.0);
			if (payoff > conti_value || (n == N))
			{
				price[i] = payoff;
				for (int k = 0; k < n; k++)
				{
					price[i] = price[i] * alpha;
				}
				break;
			}
		}
	}
	double p0 = 0;
	for (int i = 0; i < m; i++)
	{
		p0 += price[i];
	}
	p0 /= m;

	return p0;
}