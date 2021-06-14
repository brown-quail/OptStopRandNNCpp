#pragma once

#include "stdafx.h"
#include "activationFunctions.h"


double doRLSM(payoff g, double strike, 
	double*** path, int m, int d, int N,
	double discount, int K, double T,
	mt19937_64 urng);