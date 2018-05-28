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

SCIP_RETCODE Load_Original_Model(structGenCol & sGC)
{
	SCIP_PROBDATA* probdata = NULL;
	
	/* create environment */
    SCIP_CALL( SCIPcreate(&sGC.scip));
    assert(sGC.scip != NULL);
    
    /* Load plugin */
    SCIP_CALL( SCIPincludeDefaultPlugins(sGC.scip) );
    
    /** project plugins */
    //SCIP_CALL( SCIPincludePricerVRPTF (pU->scip, pU));
    
    probdata = (SCIP_PROBDATA*) &sGC;
    
    /* create empty problem */
    SCIP_CALL( SCIPcreateProb(sGC.scip, "Problem_to_Solve_for_Test", 0, 0, 0, 0, 0, 0, probdata) );
    
    // variables x_it
	vector<vector<SCIP_VAR *> > x_it;
	for(int i=0; i<sGC.d.cardJ; ++i){
		vector<SCIP_VAR *> v;
		//for(int t=sGC.d.rj[i]; t<sGC.d.dj[i]; ++t){
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_VAR * var;
			SCIPcreateVarBasic(sGC.scip, &var, ("x_it"+to_string(i)+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(sGC.scip,var);
			v.push_back(var);		
		}
		x_it.push_back(v);
	}
	sGC.varX_it = x_it;
	cout<<"x ok"<<endl;
	
    // variables z_pt et z0_pt
	vector<vector<SCIP_VAR *> > z_pt;
    vector<vector<SCIP_VAR *> > z0_pt;
	for(int p=0;p<sGC.d.nb_bp[0]; ++p){
		vector<SCIP_VAR *> v;
        vector<SCIP_VAR *> v0;
		z_pt.push_back(v);
        z0_pt.push_back(v0);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_VAR * var;
            SCIP_VAR * var0;
			z_pt[p].push_back(var);
            z0_pt[p].push_back(var0);
			SCIPcreateVarBasic(sGC.scip, &(z_pt[p][t]), ("z_pt"+to_string(p)+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
            SCIPcreateVarBasic(sGC.scip, &(z0_pt[p][t]), ("z0_pt"+to_string(p)+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(sGC.scip,z_pt[p][t]);
            SCIPaddVar(sGC.scip,z0_pt[p][t]);		
		}	
	}

	cout<<"z ok"<<endl;
    // variables y_lkt
	vector<vector<vector<SCIP_VAR *> > > y_lkt;
	for(int l=0; l<sGC.cardL; ++l){
		vector<vector<SCIP_VAR *> > v;
		y_lkt.push_back(v);
		for(int k=0; k<sGC.d.nb_bp[0]; ++k){
			vector<SCIP_VAR *> var;
			y_lkt[l].push_back(var);
			for(int t=0; t<sGC.d.cardT; ++t){
				SCIP_VAR * varr;
				y_lkt[l][k].push_back(varr);
				SCIPcreateVarBasic(sGC.scip, &y_lkt[l][k][t], ("y_lkt"+to_string(l)+to_string(k)+to_string(t)).c_str(), 0.0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_CONTINUOUS);
				SCIPaddVar(sGC.scip,y_lkt[l][k][t]);		
			}
		}	
	}
	sGC.varY_lkt = y_lkt;
    cout<<"y ok"<<endl;
    // variables y0_kt
	vector<vector<SCIP_VAR *> > y0_kt;
	for(int k=0; k<sGC.d.nb_bp[0]; ++k){
		vector<SCIP_VAR *> v;
		y0_kt.push_back(v);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_VAR * var;
			y0_kt[k].push_back(var);
			SCIPcreateVarBasic(sGC.scip, &(y0_kt[k][t]), ("y0_kt"+to_string(k)+to_string(t)).c_str(), 0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(sGC.scip,y0_kt[k][t]);		
		}	
	}
	sGC.varY0_kt = y0_kt;
	cout<<"y0 ok"<<endl;

    // variables q_t
	vector<SCIP_VAR *> q_t;
	for(int i=0; i<sGC.d.cardT; ++i){
		SCIP_VAR * var;
		q_t.push_back(var);
		SCIPcreateVarBasic(sGC.scip, &q_t[i], ("q_t"+to_string(i)).c_str(), sGC.p.qmin, sGC.p.qmax, 0.0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(sGC.scip,q_t[i]);
	}
	cout<<"q ok"<<endl;
    // fonction objectif
	for(int t=0; t<sGC.d.cardT; ++t){
		for(auto const& l : sGC.L_t[t]){
			for(auto const& k : sGC.K_l[l]){
				//cout << t <<" "<< l << " "<<k<<" "<<sGC.d.valbpt[0][k]<<endl;
				SCIPchgVarObj(sGC.scip,y_lkt[l][k][t],sGC.d.valbpt[0][k]);
				SCIPchgVarObj(sGC.scip,y0_kt[k][t],sGC.d.valbpt[0][k]);
			}
		}
	}
	cout<<"obj ok"<<endl;


	// Add constraints
    
	
	// x_it - sum_l sum_k a_il y_lkt = 0
	vector<vector<SCIP_CONS *> > cons_1;
	for(int i=0; i<sGC.d.cardJ; ++i){
		vector<SCIP_CONS *> c;
		cons_1.push_back(c);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_CONS * ct;
			cons_1[i].push_back(ct);
			SCIPcreateConsLinear(sGC.scip, &cons_1[i][t], "cons_1", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, cons_1[i][t], x_it[i][t],1);		
			for(const auto& l : sGC.L_t[t]){
				for(const auto& k : sGC.K_l[l]){
					SCIPaddCoefLinear(sGC.scip, cons_1[i][t], y_lkt[l][k][t],-sGC.a_il[i][l]);
				}
			}
			SCIPaddCons(sGC.scip, cons_1[i][t]);
		}
	}
	sGC.cons_1 = cons_1;
	cout<<"cons_1 ok"<<endl;

	// sum_t sum_l sum_k a_il y_lkt >= pi
	vector<SCIP_CONS *> cons_2;
	for(int i=0; i<sGC.d.cardJ; ++i){
		SCIP_CONS * c;
		cons_2.push_back(c);
		SCIPcreateConsLinear(sGC.scip, &cons_2[i], "cons_2", 0, 0, 0, sGC.d.pj[i], SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		for(int t=0; t<sGC.d.cardT; ++t){
			for(const auto& l : sGC.L_t[t]){
				for(const auto& k : sGC.K_l[l]){
					//cout << "a_"<<i<<l<<" : "<<sGC.a_il[i][l]<<endl;
					SCIPaddCoefLinear(sGC.scip, cons_2[i], y_lkt[l][k][t],-sGC.a_il[i][l]);
				}
			}
		}
		SCIPaddCons(sGC.scip, cons_2[i]);
	}
	sGC.cons_2 = cons_2;
	cout<<"cons_2 ok"<<endl;

	// z_pt - sum_l y_l,k1,t - sum_l y_l,k2,t = 0
	vector<vector<SCIP_CONS *> > cons_3;
	for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
		vector<SCIP_CONS *> c;
		cons_3.push_back(c);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_CONS * ct;
			cons_3[p].push_back(ct);
			SCIPcreateConsLinear(sGC.scip, &cons_3[p][t], "cons_3", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, cons_3[p][t], z_pt[p][t],1);		
			for(const auto& l : sGC.L_t[t]){
				SCIPaddCoefLinear(sGC.scip, cons_3[p][t], y_lkt[l][p*2][t],-1);
				SCIPaddCoefLinear(sGC.scip, cons_3[p][t], y_lkt[l][p*2+1][t],-1);
			}
			SCIPaddCons(sGC.scip, cons_3[p][t]);
		}
	}
	sGC.cons_3 = cons_3;
	cout<<"cons_3 ok"<<endl;

	// z0_pt - y0_k1,t - y0_k2,t = 0
	vector<vector<SCIP_CONS *> > cons_4;
	for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
		vector<SCIP_CONS *> c;
		cons_4.push_back(c);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_CONS * ct;
			cons_4[p].push_back(ct);
			SCIPcreateConsLinear(sGC.scip, &cons_4[p][t], "cons_4", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, cons_4[p][t], z0_pt[p][t],1);		
			SCIPaddCoefLinear(sGC.scip, cons_4[p][t], y0_kt[p*2][t],-1);
			SCIPaddCoefLinear(sGC.scip, cons_4[p][t], y0_kt[p*2+1][t],-1);
			SCIPaddCons(sGC.scip, cons_4[p][t]);
		}
	}
	sGC.cons_4 = cons_4;
	cout<<"cons_4 ok"<<endl;

	// sum_p z_pt + z0_pt <= 1
	vector<SCIP_CONS *> cons_5;
	for(int t=0; t<sGC.d.cardT; ++t){
		SCIP_CONS * c;
		cons_5.push_back(c);
		SCIPcreateConsLinear(sGC.scip, &cons_5[t], "cons_5", 0, 0, 0, -SCIPinfinity(sGC.scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
			SCIPaddCoefLinear(sGC.scip, cons_5[t], z_pt[p][t],1);
			SCIPaddCoefLinear(sGC.scip, cons_5[t], z0_pt[p][t],1);
		}
		SCIPaddCons(sGC.scip, cons_5[t]);
	}
	cout<<"cons_5 ok"<<endl;

	// q0 = qinit
	SCIP_CONS * cons_6;
	SCIPcreateConsLinear(sGC.scip, &cons_6, "cons_6", 0, 0, 0, sGC.p.qinit, sGC.p.qinit, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(sGC.scip, cons_6, q_t[0], 1);
	SCIPaddCons(sGC.scip, cons_6);
	cout<<"cons_6 ok"<<endl;

	// q_max - q0 >= 0
	SCIP_CONS * cons_7;
	SCIPcreateConsLinear(sGC.scip, &cons_7, "cons_7", 0, 0, 0, 0, SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(sGC.scip, cons_7, q_t[sGC.d.cardT-1], 1);
	SCIPaddCoefLinear(sGC.scip, cons_7, q_t[0], -1);
	SCIPaddCons(sGC.scip, cons_7);
	cout<<"cons_7 ok"<<endl;

	// qt+1 - qt + sum_l sum_k (gout_lk - gin_lk)*y_lkt - gin_lk*y0_kt = 0
	vector<SCIP_CONS *> cons_8;
	for(int t=0; t<sGC.d.cardT-1; ++t){
		SCIP_CONS * c;
		cons_8.push_back(c);
		SCIPcreateConsLinear(sGC.scip, &cons_8[t], "cons_8", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(sGC.scip, cons_8[t], q_t[t+1],1);
		SCIPaddCoefLinear(sGC.scip, cons_8[t], q_t[t],-1);
		for(const auto& l : sGC.L_t[t]){
			for(const auto& k : sGC.K_l[l]){
				SCIPaddCoefLinear(sGC.scip, cons_8[t], y_lkt[l][k][t],(sGC.L[l].energyDemand-sGC.d.bpt[0][k]));
				SCIPaddCoefLinear(sGC.scip, cons_8[t], y0_kt[k][t],-sGC.d.bpt[0][k]);
			}
		}
		SCIPaddCons(sGC.scip, cons_8[t]);
	}
	int t = sGC.d.cardT-1;
	SCIP_CONS * c;
	cons_8.push_back(c);
	SCIPcreateConsLinear(sGC.scip, &cons_8[t], "cons_8", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(sGC.scip, cons_8[t], q_t[t],-1);
	for(const auto& l : sGC.L_t[t]){
		for(const auto& k : sGC.K_l[l]){
			SCIPaddCoefLinear(sGC.scip, cons_8[t], y_lkt[l][k][t],(sGC.L[l].energyDemand-sGC.d.bpt[0][k]));
			SCIPaddCoefLinear(sGC.scip, cons_8[t], y0_kt[k][t],-sGC.d.bpt[0][k]);
		}
	}
	SCIPaddCons(sGC.scip, cons_8[t]);
	sGC.cons_8 = cons_8;
	cout<<"cons_8 ok"<<endl;

	// sum_l sum_k (ok - gin_lk + gout_lk)*y_lkt - sum_i bi x_it = 0
	vector<SCIP_CONS *> cons_9;
	for(int t=0; t<sGC.d.cardT; ++t){
		SCIP_CONS * c;
		cons_9.push_back(c);
		SCIPcreateConsLinear(sGC.scip, &cons_9[t], "cons_9", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		for(const auto& l : sGC.L_t[t]){
			for(const auto& k : sGC.K_l[l]){
				SCIPaddCoefLinear(sGC.scip, cons_9[t], y_lkt[l][k][t],sGC.L[l].energyDemand);
			}
		}
		for(int i=0; i<sGC.d.cardJ; ++i) SCIPaddCoefLinear(sGC.scip, cons_9[t], x_it[i][t],sGC.d.dj[i]);
		SCIPaddCons(sGC.scip, cons_9[t]);
	}
	sGC.cons_9 = cons_9;
	cout<<"cons_9 ok"<<endl;

	sGC.nbconstmodele = SCIPgetNConss(sGC.scip);
	sGC.nbvarmodele = SCIPgetNVars(sGC.scip);
    sGC.nbcolgenerated = 0;
    
	
	FILE * filed;
	filed = fopen("sol_comp", "w");
	SCIPprintOrigProblem(sGC.scip, filed, "lp", false);
	
    return SCIP_OKAY;
	
}