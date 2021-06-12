#pragma once

#include <EigenRand/Core.h>
#include "activationFunctions.h"

typedef double (*payoff)(double***, int, int, int, double);

using namespace std;
double doRLSM(payoff g, double strike, 
	double*** path, int m, int d, int N,
	double discount, int K, double T,
	mt19937_64 urng);