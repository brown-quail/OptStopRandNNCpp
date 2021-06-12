#include "payoff.h"
#include <math.h>
#include <algorithm>

double max_call(double*** x, int simNo, int n, int assetDim, double K)
{
	double payoff = 0;
	for (int i = 0; i < assetDim; i++)
	{
		payoff = max(payoff, x[simNo][i][n]);
	}
	payoff -= K;
	payoff = max(payoff, 0.0);
	return payoff;
}

double geometric_put(double*** x, int simNo, int n, int assetDim, double K)
{
	double payoff = 1;
	for (int i = 0; i < assetDim; i++)
	{
		payoff *= x[simNo][i][n];
	}
	payoff = pow(payoff, 1 / assetDim);
	payoff = K - payoff;
	payoff = max(payoff, 0.0);
	return payoff;
}

double basket_call(double*** x, int simNo, int n, int assetDim, double K)
{
	double payoff = 0;
	for (int i = 0; i < assetDim; i++)
	{
		payoff += x[simNo][i][n];
	}
	payoff /= assetDim;
	payoff -= K;
	payoff = max(payoff, 0.0);
	return payoff;
}