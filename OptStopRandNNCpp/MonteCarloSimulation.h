#pragma once
#include "stdafx.h"

double*** BlackScholes(double rf, double dividend, double vol, 
	double x0, int dimension, double T, double M, double N, mt19937_64 urng);

double*** Heston(double rf, double volvol, double ltvar, double kappa, double HestonCorr, 
	double x0, double nu0, int dimension, double T, double M, double N, mt19937_64 urng);

//double*** fractionalBlackScholes(double rf, double dividend, double vol, double H,
//	double x0, int dimension, double T, double M, double N, mt19937_64 urng);
// not implemented

void deletePath(double*** path, int m, int d);