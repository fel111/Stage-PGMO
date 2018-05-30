#ifndef struct_gencol_H
#define struct_gencol_H

#include <vector>
#include "../struct.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

using namespace std;


struct feasibleSet{
    int id;
	int releaseTime;
	int deadLine;
	float energyDemand;
	float cost;
	vector<int> tasksList;
	int timeGen;
};

struct structGenCol{
    data d;
    param p;
    SCIP *scip;
    int nbconstmodele, nbvarmodele, nbcolgenerated;
    vector<feasibleSet> L;
    int cardL;

    vector<vector<int> > L_t;

    vector<vector<int> > K_l;

    vector<vector<int> > a_il;
    // variables pour generation colonnes
    vector<vector<vector<SCIP_VAR *> > > varY_lkt;
    vector<vector<SCIP_VAR *> > varY0_kt;
    vector<vector<SCIP_VAR *> > varX_it;
    //vector<vector<SCIP_VAR *> > varZ_pt;
    //vector<vector<SCIP_VAR *> > varZ0_pt;

    //contraintes pour generation colonnes
    vector<vector<SCIP_CONS *> > cons_1;
    vector<SCIP_CONS *> cons_2;
    vector<vector<SCIP_CONS *> > cons_3;
    //vector<vector<SCIP_CONS *> > cons_4;
    vector<SCIP_CONS *> cons_8;
    vector<SCIP_CONS *> cons_9;

    SCIP_SOL * sol;

};



float cl(int k, int l, structGenCol const& sGC);



#endif
