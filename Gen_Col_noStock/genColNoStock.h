#ifndef GENCOL_NOSTOCK_H
#define GENCOL_NOSTOCK_H
#include "struct_gencol.h"
#include "../struct.h"

using namespace std;

SCIP_RETCODE Load_Original_Model(structGenCol & sGC);
float firstSol(structGenCol &sGC);
float genColNoStock(data & d,param & p,float &tps,vector<float> & demande, float &gap, string &statut);

#endif