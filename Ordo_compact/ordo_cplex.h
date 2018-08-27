#ifndef ORDO_CPLEX_H
#define ORDO_CPLEX_H

//#include <vector>
#include "../struct.h"
//using namespace std;

float ordo_cplex(data const& d,param const& p, float &tps, vector<float> &demande,string &status, float &gap);


#endif
