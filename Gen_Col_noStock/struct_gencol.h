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
	//int timeGen;
};

struct structGenCol{
    int todelete;

    float tpsSP= 0.0;
    float tpsPM= 0.0;
    float tpsPricer=0.0;

    data d;
    param p;
    SCIP *scip;
    int nbcolgenerated=0;
    int cptPricer=0;
    vector<feasibleSet> L;
    int cardL;

    vector<vector<int> > L_t;

    vector<vector<int> > K_l;

    vector<vector<int> > P_l;

    vector<vector<int> > a_il;

    vector<int> P_0;
    // variables pour generation colonnes
    vector<vector<SCIP_VAR *> > varY_lt;
    vector<vector<SCIP_VAR *> > varX_it;

    //contraintes pour generation colonnes
    vector<vector<SCIP_CONS *> > cons_1;
    vector<SCIP_CONS *> cons_2;
    vector<SCIP_CONS *> cons_3;

    //variables duales
    vector<vector<double> > w_it;
    vector<double> u_i;
    vector<vector<double> > z_it;
    vector<double> v_t;

    //int time;

    SCIP_SOL * sol;

    structGenCol(const data & dat, const param & par) : d(dat), p(par) {}
    structGenCol(){}

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
float p_t(int t, float x, structGenCol const& sGC);


#endif
