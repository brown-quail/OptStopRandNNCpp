#pragma once
#include <EigenRand/Core.h>
using namespace std;

double*** BlackScholes(double rf, double dividend, double vol, 
	double x0, int dimension, double T, double M, double N, mt19937_64 urng);

double*** Heston(double rf, double volvol, double ltvar, double kappa, double HestonCorr, 
	double x0, double nu0, int dimension, double T, double M, double N, mt19937_64 urng);

void deletePath(double*** path, int m, int d);