#ifndef MODELE_ENTIER_CPLEX_H
#define MODELE_ENTIER_CPLEX_H
#include <vector>
#include "data_struct.h"



float modele_entier_cplex(data d,param p);
float relaxation_modele_entier_cplex(data d,vector<float>& relax,param p);

#endif