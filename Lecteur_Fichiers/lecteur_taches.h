#ifndef LECTEUR_TACHES_H
#define LECTEUR_TACHES_H
//#include "data_struct.h"

using namespace std;

void lecteur_taches(string file, data &d,bool fenetre);
void generateur_taches(string file, int cardJ, int cardT, int dureemax, float consomax);
void lecteur_taches_EnergSchedInst(string file, data &d);

#endif