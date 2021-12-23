#include <iostream>
#include <vector>
#include <EigenRand/Core.h>

#include "MonteCarloSimulation.h"
#include "payoff.h"
#include "RFQI.h"
#include "RLSM.h"
#include "RRLSM.h"
#include "LSM.h"
#include <ctime>
#include <fstream>

using namespace std;


int main(void)
{
	bool isRLSM = true;
	bool isRFQI = false;
	bool isLSM = true;


	int max_iter = 30;
	double tol = 1e-4;
	double r = 0.02;
	double div = 0.0;
	double vol = 0.2;
	double x0 = 100;
	int d = 1;
	double corr = 0;
	double discount = r;

	double strike = 100;

	double T = 3;
	int N = 1000;
	int m = 20000;

	int RLSM_K = 20;
	int RFQI_K = std::max(std::min(d, 20), 5);

	int experiments = 10;

	ofstream fout;
	fout.open("OptStopRandNNCppOutput.txt", ios::trunc);
	fout << "MCS_Time";
	if (isRLSM)	fout << "\tRLSM_Time";
	if (isRFQI)	fout << "\tRFQI_Time";
	if (isLSM )	fout << "\t LSM_Time";
	if (isRLSM)	fout << "\tRLSM_price";
	if (isRFQI)	fout << "\tRFQI_price";
	if (isLSM )	fout << "\t LSM_price";
	fout << endl;

	double mcs_time, rlsm_time, rfqi_time, lsm_time;
	double RFQI_price, RLSM_price, LSM_price;

	for (int i = 0; i < experiments; i++)
	{
		cout << "[Experiment " << i << "]" << endl;
		std::random_device rd;
		int seed = rd();
		std::cout << "seed: " << seed << std::endl;
		std::mt19937_64 urng(seed);
		
		clock_t startTime = clock();
		double*** path = BlackScholes(r, div, vol, x0, d, T, m, N, urng);
		mcs_time = (double)(clock() - startTime) / CLOCKS_PER_SEC;
		cout << "Path generation completed!\telapsed(s): " << mcs_time << endl;

		if (isRFQI)
		{
			startTime = clock();
			RFQI_price = doRFQI(max_call, strike, path, m, d, N, discount, RFQI_K, T, max_iter, tol, urng);
			rfqi_time = (double)(clock() - startTime) / CLOCKS_PER_SEC;
			cout << "RFQI Price: " << RFQI_price << "\telapsed(s): " << rfqi_time << endl;
		}

		if (isRLSM)
		{
			startTime = clock();
			RLSM_price = doRLSM(max_call, strike, path, m, d, N, discount, RLSM_K, T, urng);
			rlsm_time = (double)(clock() - startTime) / CLOCKS_PER_SEC;
			cout << "RLSM Price: " << RLSM_price << "\telapsed(s): " << rlsm_time << endl;
		}

		if (isLSM)
		{
			startTime = clock();
			LSM_price = doLSM2(max_call, strike, path, m, d, N, discount, T);
			lsm_time = (double)(clock() - startTime) / CLOCKS_PER_SEC;
			cout << " LSM Price: " << LSM_price << "\telapsed(s): " << lsm_time << endl;
		}

		deletePath(path, m, d);

		fout << mcs_time;
		if (isRLSM) fout << "\t" << rlsm_time;
		if (isRFQI) fout << "\t" << rfqi_time;
		if (isLSM ) fout << "\t" << lsm_time;
		if (isRLSM) fout << "\t" << RLSM_price;
		if (isRFQI) fout << "\t" << RFQI_price;
		if (isLSM ) fout << "\t" << LSM_price;
		fout << endl;

		cout << endl;

	}

	fout.close();

	return 0;
}