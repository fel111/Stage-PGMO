#ifndef LOTSIZCONTCPX_H
#define LOTSIZCONTCPX_H

#include "struct.h"
#include <vector>

using namespace std;

float lotsizcontCPX(data const& d, param const& p, vector<float> &varia, float &tps, float &borneinf, string &status);
//float lotsizcontCPX(data const& d, float& borninf, string &status, float & tps);


#endif