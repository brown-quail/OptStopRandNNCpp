#include <iostream>
#include <vector>
#include <EigenRand/Core.h>

#include "MonteCarloSimulation.h"
#include "payoff.h"
#include "RFQI.h"
#include "RLSM.h"

using namespace std;


int main(void)
{
	int max_iter = 30;
	double tol = 1e-4;
	double r = 0.02;
	double div = 0.0;
	double vol = 0.2;
	double x0 = 100;
	int d = 5;
	double corr = 0;
	double discount = r;

	double strike = 100;

	double T = 1;
	int N = 10;
	int m = 20000;

	int RLSM_K = 20;
	int RFQI_K = std::max(std::min(d, 20), 5);

	std::random_device rd;
	int seed = rd();
	std::cout << "seed: " << seed << std::endl;
	std::mt19937_64 urng(seed);

	double*** path = BlackScholes(r, div, vol, x0, d, T, m, N, urng);
	cout << "Path generation completed!" << endl;

	//double RFQI_price = doRFQI(max_call, strike, path, m, d, N, discount, RFQI_K, T, max_iter, tol, urng);

	//cout << "RFQI Price: " << RFQI_price << endl;

	double RLSM_price = doRLSM(max_call, strike, path, m, d, N, discount, RLSM_K, T, urng);

	cout << "RLSM Price: " << RLSM_price << endl;





	deletePath(path, m, d);

	return 0;
}