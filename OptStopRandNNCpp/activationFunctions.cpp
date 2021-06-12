#include "activationFunctions.h"
#include <Eigen/Dense>

double leakyRelu(double x)
{
	return x >= 0 ? x : x * 0.1;
}

Eigen::MatrixXd leakyRelu(Eigen::MatrixXd mat)	// matrix-version overloading
{
	for (int i = 0; i < mat.rows(); i++)
		for (int j = 0; j < mat.cols(); j++)
		{
			mat(i, j) = leakyRelu(mat(i, j));
		}
	return mat;
}