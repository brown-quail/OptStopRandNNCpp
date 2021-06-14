#pragma once

#include "stdafx.h"

double doLSM2(payoff g, double strike,
	double*** path, int m, int d, int N,
	double discount, double T);