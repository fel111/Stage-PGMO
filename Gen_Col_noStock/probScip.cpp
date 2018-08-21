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

SCIP_RETCODE Load_Original_Model(structGenCol & sGC)
{
	SCIP_PROBDATA* probdata = NULL;
	
	/* create environment */
    SCIP_CALL( SCIPcreate(&sGC.scip));
    assert(sGC.scip != NULL);
    
    /* Load plugin */
    SCIP_CALL( SCIPincludeDefaultPlugins(sGC.scip) );
    if(sGC.p.aff_log_ordo == 0) SCIPsetMessagehdlr(sGC.scip,NULL);
	SCIP_CALL( SCIPchgRealParam(sGC.scip,SCIPgetParam(sGC.scip,"limits/time"),sGC.p.time_limit_ordo) );
	//SCIP_CALL( SCIPsetIntParam(sGC.scip, "presolving/maxrounds", 0));
	//SCIP_CALL( SCIPchgRealParam(sGC.scip,SCIPgetParam(sGC.scip,"numerics/epsilon"),0.0001) );
    //SCIP_CALL( SCIPchgIntParam(sGC.scip,SCIPgetParam(sGC.scip,"lp/colagelimit"),-1) ); // permet d'empecher aging (marche pas)
	//SCIP_CALL( SCIPchgLongintParam(sGC.scip,SCIPgetParam(sGC.scip,"limits/totalnodes"),1) );
	//SCIP_CALL( SCIPchgStringParam(sGC.scip,SCIPgetParam(sGC.scip,"visual/vbcfilename"),"vbcfileNoStock.vbc") );
	//SCIP_CALL( SCIPchgBoolParam(sGC.scip,SCIPgetParam(sGC.scip,"visual/dispsols"),TRUE) );


	/** project plugins */
    SCIP_CALL( includePricer(sGC) );
    
    probdata = (SCIP_PROBDATA*) &sGC;
    
    /* create empty problem */
    SCIP_CALL( SCIPcreateProb(sGC.scip, "Problem_GenCol_NoStock", 0, 0, 0, 0, 0, 0, probdata) );
    
    // variables x_it
	vector<vector<SCIP_VAR *> > x_it;
	for(int i=0; i<sGC.d.cardJ; ++i){
		vector<SCIP_VAR *> v (sGC.d.cardT,NULL);
		for(int t=sGC.d.rj[i]; t<sGC.d.dj[i]; ++t){
			SCIP_VAR * var;
			SCIPcreateVarBasic(sGC.scip, &var, ("x_it"+to_string(i)+','+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(sGC.scip,var);
			v[t] = var;	
		}
		x_it.push_back(v);
	}
	sGC.varX_it = x_it;
	
    // variables y_lt
	vector<vector<SCIP_VAR *> >  y_lt;
	for(int l=0; l<sGC.L.size(); ++l){
		vector<SCIP_VAR *> v (sGC.d.cardT, NULL);
		y_lt.push_back(v);
		for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
			SCIP_VAR * varr;
			SCIPcreateVarBasic(sGC.scip, &varr, ("y_lt"+to_string(l)+','+to_string(t)).c_str(), 0.0, SCIPinfinity(sGC.scip), sGC.L[l].cost, SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(sGC.scip,varr);
			y_lt[l][t] = varr;
		}
	}
	sGC.varY_lt = y_lt;

    

	// Add constraints
    
	
	// x_it - sum_l a_il y_lkt = 0
	vector<vector<SCIP_CONS *> > cons_1;
	for(int i=0; i<sGC.d.cardJ; ++i){
		vector<SCIP_CONS *> c;
		cons_1.push_back(c);
		for(int t=sGC.d.rj[i]; t<sGC.d.dj[i]; ++t){
			SCIP_CONS * ct;
			SCIPcreateConsLinear(sGC.scip, &ct, ("cons_1,"+to_string(i)+","+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, ct, sGC.varX_it[i][t],1);		
			for(const auto& l : sGC.L_t[t]){
				SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lt[l][t],-sGC.a_il[i][l]);
			}
			SCIP_CALL(SCIPaddCons(sGC.scip, ct));
			cons_1[i].push_back(ct);
		}
		vector<double> vd (sGC.d.dj[i]-sGC.d.rj[i],0.0);
		sGC.w_it.push_back(vd);
	}
	sGC.cons_1 = cons_1;

	// sum_t sum_l a_il y_lt >= pi
	vector<SCIP_CONS *> cons_2;
	for(int i=0; i<sGC.d.cardJ; ++i){
		SCIP_CONS * c;
		SCIPcreateConsLinear(sGC.scip, &c, ("cons_2,"+to_string(i)).c_str(), 0, 0, 0, sGC.d.pj[i], SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		
		for(int l=0; l<sGC.L.size(); ++l){
			for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
				//cout << "a_"<<i<<l<<" : "<<sGC.a_il[i][l]<<endl;
				SCIPaddCoefLinear(sGC.scip, c, sGC.varY_lt[l][t],sGC.a_il[i][l]);
			}
		}
		
		SCIPaddCons(sGC.scip, c);
		cons_2.push_back(c);
	}
	vector<double> vdj (sGC.d.cardJ,0.0);
	sGC.u_i = vdj;
	sGC.cons_2 = cons_2;
	//cout<<"cons_2 ok"<<endl;

	// sum_l -ylt >= -1
	vector<SCIP_CONS *> cons_3;
	for(int t=0; t<sGC.d.cardT; ++t){
		SCIP_CONS * ct;
		SCIPcreateConsLinear(sGC.scip, &ct, ("cons_3,"+to_string(t)).c_str(), 0, 0, 0, -1, SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  	
		for(const auto& l : sGC.L_t[t]){
			SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lt[l][t],-1);		
		}
		SCIPaddCons(sGC.scip, ct);
		cons_3.push_back(ct);
	}
	sGC.cons_3 = cons_3;
	vector<double> vdt (sGC.d.cardT,0.0);
	sGC.v_t = vdt; 
	//cout<<"cons_3 ok"<<endl;

	
    return SCIP_OKAY;
	
}