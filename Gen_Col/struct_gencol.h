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
    int todelete;
    int todelete2=0;
    double todelete3=0.0;

    float tpsSP= 0.0;
    float tpsPM= 0.0;
    float tpsPricer=0.0;

    data d;
    param p;
    SCIP *scip;
    int nbcolgenerated=0;
    int cptPricer=0;
    int cptFarkas=0;
    vector<feasibleSet> L;
    int cardL;

    vector<vector<int> > L_t;

    vector<vector<int> > K_l;

    vector<vector<int> > P_l;

    vector<vector<int> > a_il;

    vector<int> P_0;
    // variables pour generation colonnes
    vector<vector<vector<SCIP_VAR *> > > varY_lkt;
    vector<vector<SCIP_VAR *> > varY0_kt;
    vector<vector<SCIP_VAR *> > varX_it;

    //contraintes pour generation colonnes
    vector<vector<SCIP_CONS *> > cons_1;
    vector<SCIP_CONS *> cons_2;
    vector<vector<SCIP_CONS *> > cons_3;
    //vector<vector<SCIP_CONS *> > cons_4;
    vector<SCIP_CONS *> cons_8;
    vector<SCIP_CONS *> cons_9;

    //variables duales
    vector<vector<double> > alpha_it;
    vector<double> beta_i;
    vector<vector<double> > gamma_pt;
    vector<double> delta_t;
    vector<double> phi_t;

    //farkas
    vector<vector<double> > farkas1_it;
    vector<double> farkas2_i;
    vector<vector<double> > farkas3_pt;
    vector<double> farkas8_t;
    vector<double> farkas9_t;

    //int time;

    SCIP_SOL * sol;

};



float cl(int k, int l, structGenCol const& sGC);
void addP_0(structGenCol & sGC);
void addL_t(feasibleSet const& l, structGenCol & sGC);
void addA_il(feasibleSet const& l, structGenCol & sGC);
void addSetK_l(feasibleSet const& l, structGenCol & sGC);
void addP_l(feasibleSet const& l, structGenCol & sGC);
int checkSet(feasibleSet const& l, structGenCol const& sGC);
void affAllSet(structGenCol const& sGC);
void affL_t(structGenCol const& sGC);
void affK_l(structGenCol const& sGC);
void modifPwlCplex(structGenCol & sGC);


#endif
