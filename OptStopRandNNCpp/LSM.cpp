#include "LSM.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

double doLSM2(payoff g, double strike,
	double*** path, int m, int d, int N,
	double discount, double T)
{
	// 여기서는 형평성을 맞추기 위해 m+1 ~ 2m까지의 path만을 사용한다.
	// 2차식 polynomial regression 만 한다.

	double alpha = exp(-discount * T / N);

	int* stopping = new int[m];
	for (int i = 0; i < m; i++)
	{
		stopping[i] = N;
	}

	for (int n = N - 1; n > 0; n--)
	{
		Eigen::MatrixXd A(m, 1 + d + d * (d + 1) / 2);
		Eigen::MatrixXd theta(1 + d + d * (d + 1) / 2, 1);
		Eigen::MatrixXd b(m, 1);
		for (int i = 0; i < m; i++)
		{
			double cont = g(path, i + m, stopping[i], d, strike);
			for (int k = 0; k < stopping[i] - n; k++)
			{
				cont *= alpha;
			}
			b(i, 0) = cont;

			int nowIdx = 0;
			A(i, nowIdx++) = 1;
			for (int j = 0; j < d; j++)
			{
				A(i, nowIdx++) = path[i][j][n];
			}
			for (int j = 0; j < d; j++)
			{
				for (int k = j; k < d; k++)
				{
					A(i, nowIdx++) = path[i][j][n] * path[i][k][n];
				}
			}
			// 1 + d + (d + d-1 + d-2 + ... + 1) = 1 + d + d*(d+1)/2
		}

		theta = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
		Eigen::MatrixXd contValue = A * theta;
		for (int i = 0; i < m; i++)
		{
			double stopValue = g(path, i + m, n, d, strike);
			if (stopValue >= contValue(i))
			{
				stopping[i] = n;
			}
		}
	}

	// when n = 0
	double contValue = 0;
	double stopValue = g(path, m, 0, d, strike);
	for (int i = 0; i < m; i++)
	{
		double cont = g(path, i + m, stopping[i], d, strike);
		for (int k = 0; k < stopping[i]; k++)
		{
			cont *= alpha;
		}
		contValue += cont;
	}
	contValue /= m;

	double p0 = max(stopValue, contValue);

	delete[] stopping;

	return p0;
}