#ifndef struct_pricer_H
#define struct_pricer_H

#include <vector>
#include "../struct.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <ilcplex/ilocplex.h>

using namespace std;

bool verifSol(structGenCol const& sGC);
double verifCoutReduit(structGenCol const& sGC, IloNumArray ui, int t);
double verifCoutReduit(structGenCol const& sGC, vector<int> ui, int t);
SCIP_RETCODE addObjectColumnInModel (structGenCol &sGC,IloNumArray const& valUi,int timedd, double objvalue);
SCIP_RESULT Pr_SP1(structGenCol &sGC);
static SCIP_DECL_PRICERREDCOST(pricerRedcost);
static SCIP_DECL_PRICERINITSOL(pricerInitsolSP);
SCIP_RETCODE includePricer(structGenCol &sGC);







#endif