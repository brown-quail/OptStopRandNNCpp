#pragma once
#include "stdafx.h"

double leakyRelu(double x);

Eigen::MatrixXd leakyRelu(Eigen::MatrixXd mat);	// matrix-version overloading