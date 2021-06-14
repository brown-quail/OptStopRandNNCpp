#include "RLSM.h"
#include <Eigen/Dense>
#include <math.h>
#include <iostream>
#include <Eigen/SVD>

double doRRLSM(payoff g, double strike, double*** path, int m, int d, int N,
	double discount, int K, double T,
	mt19937_64 urng)
{
	Eigen::MatrixXd Ax = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, d, urng);
	Eigen::MatrixXd Ah = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, K - 1, urng);
	Eigen::MatrixXd b = Eigen::Rand::normal<Eigen::MatrixXd>(K - 1, 1, urng);

	Eigen::MatrixXd theta(K, 1);

	theta.setZero();

	double** pni = new double* [N + 1];
	for (int i = 0; i < N + 1; i++)
	{
		pni[i] = new double[2 * m];
	}

	Eigen::MatrixXd phi(K, 1);
	Eigen::MatrixXd x(d, 1);

	Eigen::MatrixXd* phiM = new Eigen::MatrixXd[N];
	for (int n = 0; n < N; n++)
	{
		phiM[n].resize(K, 2 * m);
		phiM[n].setConstant(0);
	}
	Eigen::MatrixXd pniM(2 * m, 1);

	double alpha = exp(-discount * T / N);

	for (int i = 0; i < 2 * m; i++)
	{
		pni[N][i] = g(path, i, N, d, strike);
	}

	for (int i = 0; i < 2 * m; i++)
	{
		Eigen::MatrixXd h(K - 1, 1);
		h.setConstant(0);	// h(0) = 0
		for (int n = 1; n < N; n++)
		{
			for (int j = 0; j < d; j++)
			{
				x(j, 0) = path[i][j][n];
			}
			h = leakyRelu(Ax * x + Ah * h + b);
			Eigen::MatrixXd temp(K - 1, 1);
			temp = h;
			temp.conservativeResize(K, 1);
			temp(K - 1, 0) = 1;
			phi = temp;
			phiM[n].col(i) = temp;
		}
	}

	for (int n = N - 1; n > 0; n--)
	{
		for (int i = 0; i < 2 * m; i++)
		{

			if (1)// (i < m)
			{
				pniM(i) = pni[n + 1][i];
			}
		}
		Eigen::MatrixXd phiMT = phiM[n].leftCols(m).transpose();	// m x K
		// solve phiMT * theta = alpha * pniM
		// dimension (m x K) * (K, 1) = (m, 1)
		theta = phiMT.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(alpha * pniM.topRows(m));
		//theta = phiMT.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(alpha * pniM.topRows(m));
		for (int i = 0; i < 2 * m; i++)
		{
			pni[n][i] = g(path, i, n, d, strike) >= (theta.transpose() * phiM[n].col(i))(0, 0) ? g(path, i, n, d, strike) : alpha * pni[n + 1][i];
		}
	}

	double p0 = 0;
	for (int i = 0; i < m; i++)
	{
		p0 += pni[1][i + m];
	}
	p0 /= m;

	p0 = max(g(path, 0, 0, d, strike), p0);

	delete[] phiM;

	for (int i = 0; i < N + 1; i++)
	{
		delete[] pni[i];
	}
	delete[] pni;

	return p0;
}