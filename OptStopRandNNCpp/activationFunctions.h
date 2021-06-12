#pragma once
#include <Eigen/Dense>

double leakyRelu(double x);

Eigen::MatrixXd leakyRelu(Eigen::MatrixXd mat);	// matrix-version overloading