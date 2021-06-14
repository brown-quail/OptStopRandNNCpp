#pragma once

#include <Eigen/Dense>
#include <EigenRand/Core.h>

typedef double (*payoff)(double***, int, int, int, double);

using namespace std;