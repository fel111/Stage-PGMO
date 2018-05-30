#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
//#include <limits>
//#include <math.h>
//#include <algorithm>
//#include <chrono>
//#include <deque>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "../struct.h"
#include "struct_gencol.h"
#include "pricer.h"

using namespace std;

bool verifSol(structGenCol const& sGC){
    for(int j=0; j<sGC.d.cardJ; ++j){
        int debut=0;
        double duree=0.0;
        int fin;
        while(SCIPisEQ(sGC.scip,SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][debut]),0)){ ++debut; }
        if(debut<sGC.d.rj[j]){ cout<<"tache "<<j<<" debut "<<debut<<"<"<<sGC.d.rj[j]<<endl; return false; }
        for(int k=debut; k<sGC.d.cardT; ++k){
            if(SCIPisEQ(sGC.scip,SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][k]),1)) fin = k;
            duree += SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][k]);
            //cout << "duree "<<duree<<endl;
        }
        if(fin>=sGC.d.dj[j]){ cout<<"fin depassée pour tache "<<j<<endl; return false; }
        if(SCIPisLT(sGC.scip,duree,sGC.d.pj[j])){ cout<<"tache "<<j<<"  duree "<<duree<<"<"<<sGC.d.pj[j]<<endl; return false; }
    }
    return true;
}




/*! \fn SCIP_DECL_PRICERINITSOL(pricerInitsolSP1)
    \brief solving process initialization method of variable pricer (called when branch and bound process is about to begin) here the pointers to the transformed variables and constraints are collected to ensure that the later procedures will use the right links
 */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolSP)
{
	structGenCol* pbdata; // the data of the pricer
	SCIP_PROBDATA* probdata;

	assert(scip != NULL);
	assert(pricer != NULL);

	probdata = SCIPgetProbData(scip);
	pbdata = (structGenCol *) probdata;
	// get transformed variables
    for(int t=0; t<pbdata->d.cardT; ++t){
		for(const auto& l : pbdata->L_t[t]){
			for(const auto& k : pbdata->K_l[l]){
	            SCIPgetTransformedVar(scip, pbdata->varY_lkt[l][k][t],&(pbdata->varY_lkt[l][k][t]);
			    assert(pbdata->varY_lkt[l][k][t],&(pbdata->varY_lkt[l][k][t]) != NULL);
            }
		}
	}

	for(int i = 0 ; i < pbdata->instance.getNbOfTasks() ; i++){
		for(int j = pbdata->instance.getTasksList().at(i).getReleaseTime() ; j < pbdata->instance.getTasksList().at(i).getDeadLine() ;j++){
			SCIPgetTransformedVar(scip, pbdata->tabVarXit.at(i).at(j-pbdata->instance.getTasksList().at(i).getReleaseTime()),
					&(pbdata->tabVarXit.at(i).at(j-pbdata->instance.getTasksList().at(i).getReleaseTime())));
			assert(pbdata->tabVarXit.at(i).at(j-pbdata->instance.getTasksList().at(i).getReleaseTime()) != NULL);
		}
	}
	// get transformed constraints
	int duration = pbdata->instance.getDeadLine()-pbdata->instance.getReleaseTime();
	int fois = 0;

	for(int i = 0 ; i < pbdata->instance.getNbOfTasks(); i++){
		for(int j = 0 ; j < duration; j++){
			if(fois==0){
				SCIPgetTransformedCons(scip,pbdata->tabConsNbSets[j],&(pbdata->tabConsNbSets[j]));
				assert(pbdata->tabConsNbSets[j] != NULL);
			}
			if(j>=pbdata->instance.getTasksList().at(i).getReleaseTime()&&j<pbdata->instance.getTasksList().at(i).getDeadLine()){
				SCIPgetTransformedCons(scip,pbdata->tabConsLinkVars[i][j],&(pbdata->tabConsLinkVars[i][j]));
				assert(pbdata->tabConsLinkVars[i][j] != NULL);
			}
		}
		SCIPgetTransformedCons(scip,pbdata->tabConsTasksDuration[i],&(pbdata->tabConsTasksDuration[i]));
		assert(pbdata->tabConsTasksDuration[i] != NULL);
		fois++;
	}
	return SCIP_OKAY;
}




/*! \brief Method creating the pricer, activates it and includes it in SCIP
 *
 \param[in,out] structGen Reference to structGenCol the entire problem data structure
 *
 * \return status: a variable of type SCIP_RETCODE, equal to SCIP_OKAY if everything went well
 */
SCIP_RETCODE includePricer(structGenCol &sGC){
	SCIP_PRICER* pricer;
	SCIP_PRICERDATA* pricerdata;
	pricerdata = (SCIP_PRICERDATA*) &sGC;
	SCIP_RETCODE status;
	pricer = NULL;
	assert (pricerdata != NULL);
	SCIP_CALL(status =  SCIPincludePricerBasic(sGC.scip, &pricer, PRICER_NAMESP1, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY, pricerRedcost, 0, pricerdata));
	assert(pricer != NULL);
	SCIP_CALL(status = SCIPsetPricerInitsol(sGC.scip, pricer, pricerInitsolSP));
        return status;
}