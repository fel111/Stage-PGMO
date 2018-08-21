#ifndef GENCOL_NOSTOCK_H
#define GENCOL_NOSTOCK_H
#include "struct_gencol.h"
#include "../struct.h"

using namespace std;

SCIP_RETCODE Load_Original_Model(structGenCol & sGC);
float firstSol(structGenCol &sGC);
float genColNoStock(const data &d, const param &p,float &tps, vector<float> & demande);

#endif