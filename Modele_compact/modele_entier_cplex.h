#ifndef MODELE_ENTIER_CPLEX_H
#define MODELE_ENTIER_CPLEX_H
#include <vector>
#include "data_struct.h"



float modele_entier_cplex(data const& d, param const& p, float &tps, float &borneinf, string &status);
float relaxation_modele_entier_cplex(data const& d,vector<float>& relax, param const& p, float &tps, float &borneinf, string &status);

#endif