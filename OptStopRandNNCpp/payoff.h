#pragma once

#include "stdafx.h"

// defines max call payoff when stopped.
double max_call(double*** x, int simNo, int n, int assetDim, double K = 100.0);
double geometric_put(double*** x, int simNo, int n, int assetDim, double K = 100.0);
double basket_call(double*** x, int simNo, int n, int assetDim, double K = 100.0);