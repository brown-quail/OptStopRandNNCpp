#include "MonteCarloSimulation.h"
#include <iostream>
#include <Eigen/Dense>
using namespace std;

bool useAntithetic = true;

// generate scenario path (number of paths(2*M) x path dimension(dimension) x time(N))
double*** BlackScholes(double rf, double dividend, double vol, 
	double x0, int dimension, double T, double M, double N, mt19937_64 urng)
{
	/*
	Eigen::MatrixXd corrMat(dimension, dimension);
	corrMat.setConstant(corr);	// make sure: - 1.0 / (dimension - 1) <= corr <= 1
	for (int i = 0; i < dimension; i++)
	{
		corrMat(i, i) = 1;
	}
	Eigen::MatrixXd cho(corrMat.llt().matrixL());	// Cholesky Decomposition
	*/

	double dt = T / N;

	// vector allocation and set intial value
	double*** path = new double** [2 * M];
	for (int i = 0; i < 2 * M; i++)
	{
		path[i] = new double* [dimension];
		for (int j = 0; j < dimension; j++)
		{
			path[i][j] = new double[N + 1];
			path[i][j][0] = x0;
		}
	}

	// stock price simulation
	// use antithetic variates
	for (int i = 0; i < 2 * M; i++)
	{
		Eigen::MatrixXd mat = Eigen::Rand::normal<Eigen::MatrixXd>(dimension, N, urng);
		Eigen::MatrixXd epsilon = mat;
		// Eigen::MatrixXd epsilon = cho * mat;
		for (int j = 0; j < dimension; j++)
		{
			for (int k = 1; k < N + 1; k++)
			{
				path[i][j][k] = path[i][j][k - 1] * exp((rf - dividend - vol * vol / 2.0) * dt + vol * sqrt(dt) * epsilon(j, k - 1));
				if (useAntithetic)
				{
					path[i + 1][j][k] = path[i + 1][j][k - 1] * exp((rf - dividend - vol * vol / 2.0) * dt - vol * sqrt(dt) * epsilon(j, k - 1));
				}
			}
		}
		if (useAntithetic)
		{
			i++;
		}
	}

	return path;
}

double*** Heston(double rf, double volvol, double ltvar, double kappa, double HestonCorr,
	double x0, double nu0, int dimension, double T, double M, double N, mt19937_64 urng)
{
	double dt = T / N;

	// vector allocation and set intial value
	double*** path = new double** [2 * M];
	double*** v = new double** [2 * M];
	for (int i = 0; i < 2 * M; i++)
	{
		path[i] = new double* [dimension];
		v[i] = new double* [dimension];
		for (int j = 0; j < dimension; j++)
		{
			path[i][j] = new double[N + 1];
			v[i][j] = new double[N + 1];
			path[i][j][0] = x0;
			v[i][j][0] = nu0;
		}
	}

	Eigen::MatrixXd corrMat(2, 2);
	corrMat(0, 0) = 1; corrMat(1, 1) = 1;
	corrMat(0, 1) = HestonCorr; corrMat(1, 0) = HestonCorr;
	Eigen::MatrixXd cho(corrMat.llt().matrixL());


	// stock price simulation
	// use antithetic variates
	for (int i = 0; i < 2 * M; i++)
	{
		// Eigen::MatrixXd epsilon = cho * mat;
		for (int j = 0; j < dimension; j++)
		{
			Eigen::MatrixXd mat = Eigen::Rand::normal<Eigen::MatrixXd>(2, N, urng);
			Eigen::MatrixXd epsilon = cho * mat;
			for (int k = 1; k < N + 1; k++)
			{
				v[i][j][k] = v[i][j][k - 1] - kappa * (v[i][j][k - 1] - ltvar) * dt 
					+ volvol * sqrt(v[i][j][k - 1]) * epsilon(0, k - 1);
				path[i][j][k] = path[i][j][k - 1] + rf * path[i][j][k - 1] * dt 
					+ sqrt(v[i][j][k - 1]) * path[i][j][k - 1] * epsilon(1, k - 1);
				if (useAntithetic)
				{
					v[i + 1][j][k] = v[i + 1][j][k - 1] - kappa * (v[i + 1][j][k - 1] - ltvar) * dt 
						+ volvol * sqrt(v[i + 1][j][k - 1]) * epsilon(0, k - 1);
					path[i + 1][j][k] = path[i + 1][j][k - 1] + rf * path[i + 1][j][k - 1] * dt
						+ sqrt(v[i + 1][j][k - 1]) * path[i + 1][j][k - 1] * epsilon(1, k - 1);
				}
			}
		}
		if (useAntithetic)
		{
			i++;
		}
	}

	deletePath(v, M, dimension);

	return path;
}


void deletePath(double*** path, int m, int d)
{
	for (int i = 0; i < 2 * m; i++)
	{
		for (int j = 0; j < d; j++)
		{
			delete[] path[i][j];
		}
		delete[] path[i];
	}
	delete[] path;

	return;
}