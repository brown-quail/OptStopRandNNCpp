#pragma once

#include "stdafx.h"
#include "activationFunctions.h"

double doRFQI(payoff g, double strike, 
	double*** path, int m, int d, int N,
	double discount, int K, double T,
	int max_iter, double tol, mt19937_64 urng);