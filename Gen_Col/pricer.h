#ifndef struct_pricer_H
#define struct_pricer_H

#include <vector>
#include "../struct.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <ilcplex/ilocplex.h>

using namespace std;

bool verifSol(structGenCol const& sGC);

SCIP_RETCODE addObjectColumnInModel (structGenCol &sGC,IloNumArray valUi,IloNumArray valVk,int timedd);
SCIP_RESULT Pr_SP1(structGenCol &sGC);
SCIP_RESULT Pr_farkas(structGenCol &sGC);
static SCIP_DECL_PRICERREDCOST(pricerRedcost);
static SCIP_DECL_PRICERINITSOL(pricerInitsolSP);
SCIP_RETCODE includePricer(structGenCol &sGC);
double verifCoutReduit(structGenCol const& sGC, vector<int> ui, int k, int t);
double verifCoutReduit(structGenCol const& sGC, IloNumArray ui, int k, int t);







#endif