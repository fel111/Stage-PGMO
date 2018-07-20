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
    SCIPsetMessagehdlr(sGC.scip,NULL);
	//SCIP_CALL( SCIPchgRealParam(sGC.scip,SCIPgetParam(sGC.scip,"numerics/epsilon"),0.0001) );
    //SCIP_CALL( SCIPchgIntParam(sGC.scip,SCIPgetParam(sGC.scip,"lp/colagelimit"),-1) ); // permet d'empecher aging (marche pas)
	//SCIP_CALL( SCIPchgLongintParam(sGC.scip,SCIPgetParam(sGC.scip,"limits/totalnodes"),1) );
	//SCIP_CALL( SCIPchgStringParam(sGC.scip,SCIPgetParam(sGC.scip,"visual/vbcfilename"),"vbcfile.vbc") );
	//SCIP_CALL( SCIPchgBoolParam(sGC.scip,SCIPgetParam(sGC.scip,"visual/dispsols"),TRUE) );


	/** project plugins */
    SCIP_CALL( includePricer(sGC) );
    
    probdata = (SCIP_PROBDATA*) &sGC;
    
    /* create empty problem */
    SCIP_CALL( SCIPcreateProb(sGC.scip, "Problem_to_Solve_for_Test", 0, 0, 0, 0, 0, 0, probdata) );
    
    // variables x_it
	vector<vector<SCIP_VAR *> > x_it;
	for(int i=0; i<sGC.d.cardJ; ++i){
		vector<SCIP_VAR *> v (sGC.d.cardT,NULL);
		for(int t=sGC.d.rj[i]; t<sGC.d.dj[i]; ++t){
			SCIP_VAR * var;
			SCIPcreateVarBasic(sGC.scip, &var, ("x_it"+to_string(i)+','+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			//SCIPchgVarBranchPriority(sGC.scip,var,2); //beaucoup trop long!!! (10 noeuds -> 300 noeuds)
			SCIPaddVar(sGC.scip,var);
			v[t] = var;	
		}
		x_it.push_back(v);
	}
	sGC.varX_it = x_it;
	//cout<<"x ok"<<endl;
	

    // variables z_pt et z0_pt
	vector<vector<SCIP_VAR *> > z_pt;
    vector<vector<SCIP_VAR *> > z0_pt;
	for(int p=0;p<(sGC.d.nb_bp[0]/2); ++p){
		vector<SCIP_VAR *> v;
        vector<SCIP_VAR *> v0;
		z_pt.push_back(v);
        z0_pt.push_back(v0);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_VAR * var;
            SCIP_VAR * var0;
			z_pt[p].push_back(var);
            z0_pt[p].push_back(var0);
			SCIPcreateVarBasic(sGC.scip, &(z_pt[p][t]), ("z_pt"+to_string(p)+','+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
            SCIPcreateVarBasic(sGC.scip, &(z0_pt[p][t]), ("z0_pt"+to_string(p)+','+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(sGC.scip,z_pt[p][t]);
            SCIPaddVar(sGC.scip,z0_pt[p][t]);		
		}	
	}
	//cout<<"z ok"<<endl;


    // variables y_lkt
	vector<vector<vector<SCIP_VAR *> > > y_lkt;
	for(int l=0; l<sGC.L.size(); ++l){
		//vector<vector<SCIP_VAR *> > v;
		vector<vector<SCIP_VAR *> > v (sGC.d.nb_bp[0], vector<SCIP_VAR *> (sGC.d.cardT, NULL));
		y_lkt.push_back(v);
		for(const auto& k : sGC.K_l[l]){
			//vector<SCIP_VAR *> var (sGC.d.cardT, NULL);
			//y_lkt[l].push_back(var);
			//cout << "l releasetime " << sGC.L[l].releaseTime << endl;
			for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
				SCIP_VAR * varr;
				SCIPcreateVarBasic(sGC.scip, &varr, ("y_lkt"+to_string(l)+','+to_string(k)+','+to_string(t)).c_str(), 0.0, SCIPinfinity(sGC.scip), sGC.d.valbpt[0][k], SCIP_VARTYPE_CONTINUOUS);
				SCIPaddVar(sGC.scip,varr);
				//cout << "l, k, t : " << l << k << t << endl;
				y_lkt[l][k][t] = varr;
			}
		}	
	}
	sGC.varY_lkt = y_lkt;
    //cout<<"y ok"<<endl;

    // variables y0_kt
	/*vector<vector<SCIP_VAR *> > y0_kt;
	for(int k=0; k<sGC.d.nb_bp[0]; ++k){
		//vector<SCIP_VAR *> v (sGC.d.cardT, NULL);
		vector<SCIP_VAR *> v;
		y0_kt.push_back(v);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_VAR * var;
			SCIPcreateVarBasic(sGC.scip, &var, ("y0_kt"+to_string(k)+','+to_string(t)).c_str(), 0, SCIPinfinity(sGC.scip), sGC.d.valbpt[0][k], SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(sGC.scip,var);	
			y0_kt[k].push_back(var);	
		}	
	}
	sGC.varY0_kt = y0_kt;*/
	//cout<<"y0 ok"<<endl;

	// variables y0_kt ajoutee seulement sur les breakpoints "admissibles"
	vector<vector<SCIP_VAR *> > y0_kt;
	for(int k=0; k<sGC.d.nb_bp[0]; ++k){
		if(sGC.P_0[k/2] == 1){
			vector<SCIP_VAR *> v1;
			y0_kt.push_back(v1);
			for(int t=0; t<sGC.d.cardT; ++t){
				SCIP_VAR * var1;
				SCIPcreateVarBasic(sGC.scip, &var1, ("y0_kt"+to_string(k)+','+to_string(t)).c_str(), 0, SCIPinfinity(sGC.scip), sGC.d.valbpt[0][k], SCIP_VARTYPE_CONTINUOUS);
				SCIPaddVar(sGC.scip,var1);
				y0_kt[k].push_back(var1);
			}
		}
		else{
			vector<SCIP_VAR *> v3 (sGC.d.cardT, NULL);
			y0_kt.push_back(v3);
		}
	}
	sGC.varY0_kt = y0_kt;

    // variables q_t
	vector<SCIP_VAR *> q_t;
	for(int i=0; i<=sGC.d.cardT; ++i){
		SCIP_VAR * var;
		q_t.push_back(var);
		SCIPcreateVarBasic(sGC.scip, &q_t[i], ("q_t"+to_string(i)).c_str(), sGC.p.qmin, sGC.p.qmax, 0.0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(sGC.scip,q_t[i]);
	}
	//cout<<"q ok"<<endl;


	/*for(int l=0; l<sGC.L.size(); ++l){
		for(const auto& k : sGC.K_l[l]){
			for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
			
				//cout << "Y"<<l<<k<<t<<" : "<<sGC.varY_lkt[l][k][t] << endl;
			}
		}
	}*/


	// Add constraints
    
	
	// x_it - sum_l sum_k a_il y_lkt = 0
	vector<vector<SCIP_CONS *> > cons_1;
	for(int i=0; i<sGC.d.cardJ; ++i){
		vector<SCIP_CONS *> c;
		cons_1.push_back(c);
		//for(int t=0; t<sGC.d.cardT; ++t){
		for(int t=sGC.d.rj[i]; t<sGC.d.dj[i]; ++t){
			SCIP_CONS * ct;
			SCIPcreateConsLinear(sGC.scip, &ct, ("cons_1,"+to_string(i)+","+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, ct, sGC.varX_it[i][t],1);		
			for(const auto& l : sGC.L_t[t]){
				for(const auto& k : sGC.K_l[l]){
					//if((sGC.L[l].releaseTime<=t)&&(t<sGC.L[l].deadLine)){
						//cout << "i, t, rj, dj : "<<i<<t<<sGC.d.rj[i]<<sGC.d.dj[i]<<endl;
						//cout << "varY_lk : " <<sGC.varY_lkt[l][k][t] << endl;
						SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lkt[l][k][t],-sGC.a_il[i][l]);
					//}
				}
			}
			SCIP_CALL(SCIPaddCons(sGC.scip, ct));
			cons_1[i].push_back(ct);
		}
		vector<double> vd (sGC.d.dj[i]-sGC.d.rj[i],0.0);
		sGC.alpha_it.push_back(vd);
	}
	sGC.cons_1 = cons_1;
	sGC.farkas1_it = sGC.alpha_it;
	//cout<<"cons_1 ok"<<endl;

	// sum_t sum_l sum_k a_il y_lkt >= pi
	vector<SCIP_CONS *> cons_2;
	for(int i=0; i<sGC.d.cardJ; ++i){
		SCIP_CONS * c;
		SCIPcreateConsLinear(sGC.scip, &c, ("cons_2,"+to_string(i)).c_str(), 0, 0, 0, sGC.d.pj[i], SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		for(int t=sGC.d.rj[i]; t<sGC.d.dj[i]; ++t){
			for(const auto& l : sGC.L_t[t]){
				for(const auto& k : sGC.K_l[l]){
					//cout << "a_"<<i<<l<<" : "<<sGC.a_il[i][l]<<endl;
					SCIPaddCoefLinear(sGC.scip, c, sGC.varY_lkt[l][k][t],sGC.a_il[i][l]);
				}
			}
		}
		SCIPaddCons(sGC.scip, c);
		cons_2.push_back(c);
	}
	vector<double> vdj (sGC.d.cardJ,0.0);
	sGC.beta_i = vdj;
	sGC.farkas2_i = vdj;
	sGC.cons_2 = cons_2;
	//cout<<"cons_2 ok"<<endl;

	// z_pt - sum_l y_l,k1,t - sum_l y_l,k2,t = 0
	vector<vector<SCIP_CONS *> > cons_3;
	for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
		vector<SCIP_CONS *> c;
		cons_3.push_back(c);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_CONS * ct;
			SCIPcreateConsLinear(sGC.scip, &ct, ("cons_3,"+to_string(p)+","+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, ct, z_pt[p][t],1);		
			for(const auto& l : sGC.L_t[t]){
				if(sGC.P_l[l][p] == 1){
					SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lkt[l][p*2][t],-1);
					SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lkt[l][p*2+1][t],-1);
				}
			}
			SCIPaddCons(sGC.scip, ct);
			cons_3[p].push_back(ct);
		}
		vector<double> vd (sGC.d.cardT,0.0);
		sGC.gamma_pt.push_back(vd);
	}
	sGC.cons_3 = cons_3;
	sGC.farkas3_pt = sGC.gamma_pt;
	//cout<<"cons_3 ok"<<endl;

	// z0_pt - y0_k1,t - y0_k2,t = 0
	vector<vector<SCIP_CONS *> > cons_4;
	for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
		vector<SCIP_CONS *> c;
		cons_4.push_back(c);
		for(int t=0; t<sGC.d.cardT; ++t){
			SCIP_CONS * ct;
			SCIPcreateConsLinear(sGC.scip, &ct, ("cons_4,"+to_string(p)+","+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, ct, z0_pt[p][t],1);
			if(sGC.P_0[p]==1){	
				SCIPaddCoefLinear(sGC.scip, ct, sGC.varY0_kt[p*2][t],-1);
				SCIPaddCoefLinear(sGC.scip, ct, sGC.varY0_kt[p*2+1][t],-1);
			}
			SCIPaddCons(sGC.scip, ct);
			cons_4[p].push_back(ct);
		}
	}
	sGC.cons_4 = cons_4;
	//cout<<"cons_4 ok"<<endl;

	// sum_p z_pt + z0_pt <= 1
	vector<SCIP_CONS *> cons_5;
	for(int t=0; t<sGC.d.cardT; ++t){
		SCIP_CONS * c;
		cons_5.push_back(c);
		SCIPcreateConsLinear(sGC.scip, &cons_5[t], ("cons_5,"+to_string(t)).c_str(), 0, 0, 0, -SCIPinfinity(sGC.scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
			SCIPaddCoefLinear(sGC.scip, cons_5[t], z_pt[p][t],1);
			SCIPaddCoefLinear(sGC.scip, cons_5[t], z0_pt[p][t],1);
		}
		SCIPaddCons(sGC.scip, cons_5[t]);
	}
	//cout<<"cons_5 ok"<<endl;

	// q0 = qinit
	SCIP_CONS * cons_6;
	SCIPcreateConsLinear(sGC.scip, &cons_6, "cons_6", 0, 0, 0, sGC.p.qinit, sGC.p.qinit, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(sGC.scip, cons_6, q_t[0], 1);
	SCIPaddCons(sGC.scip, cons_6);
	//cout<<"cons_6 ok"<<endl;

	// q_max - q0 >= 0
	SCIP_CONS * cons_7;
	SCIPcreateConsLinear(sGC.scip, &cons_7, "cons_7", 0, 0, 0, 0, SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(sGC.scip, cons_7, q_t[sGC.d.cardT], 1);
	SCIPaddCoefLinear(sGC.scip, cons_7, q_t[0], -1);
	SCIPaddCons(sGC.scip, cons_7);
	//cout<<"cons_7 ok"<<endl;

	// qt+1 - qt + sum_l sum_k (gout_lk - gin_lk)*y_lkt - gin_lk*y0_kt = 0
	vector<SCIP_CONS *> cons_8;
	for(int t=0; t<sGC.d.cardT; ++t){
		SCIP_CONS * c;
		SCIPcreateConsLinear(sGC.scip, &c, ("cons_8,"+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(sGC.scip, c, q_t[t+1],1);
		SCIPaddCoefLinear(sGC.scip, c, q_t[t],-1);
		for(const auto& l : sGC.L_t[t]){
			for(const auto& k : sGC.K_l[l]){
				//cout<<"contrainte 8 : sGC.L[l].energyDemand = "<< sGC.L[l].energyDemand  <<"  sGC.d.bpt[0][k] = "<<sGC.d.bpt[0][k]<<endl;
				SCIPaddCoefLinear(sGC.scip, c, sGC.varY_lkt[l][k][t],(sGC.L[l].energyDemand-sGC.d.bpt[0][k]));
				//if(sGC.L[l].energyDemand<sGC.d.bpt[0][k]) SCIPaddCoefLinear(sGC.scip, c, y0_kt[k][t],sGC.L[l].energyDemand-sGC.d.bpt[0][k]);
				//if(sGC.L[l].tasksList.size() == 0) SCIPaddCoefLinear(sGC.scip, c, y0_kt[k][t],-sGC.d.bpt[0][k]);
				//SCIPaddCoefLinear(sGC.scip, c, y0_kt[k][t],-sGC.d.bpt[0][k]);
			}
		}
		for(int k=0; k<sGC.d.nb_bp[0];++k){
			if(sGC.P_0[k/2]==1) SCIPaddCoefLinear(sGC.scip, c, y0_kt[k][t],-sGC.d.bpt[0][k]);
		}
		SCIPaddCons(sGC.scip, c);
		cons_8.push_back(c);
	}
	vector<double> vdt (sGC.d.cardT,0.0);
	sGC.delta_t = vdt;
	sGC.farkas8_t = vdt;
	/*int t = sGC.d.cardT-1;
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
	SCIPaddCons(sGC.scip, cons_8[t]);*/
	sGC.cons_8 = cons_8;
	//cout<<"cons_8 ok"<<endl;

	// sum_l sum_k (ok - gin_lk + gout_lk)*y_lkt - sum_i bi x_it = 0
	vector<SCIP_CONS *> cons_9;
	for(int t=0; t<sGC.d.cardT; ++t){
		SCIP_CONS * c;
		SCIPcreateConsLinear(sGC.scip, &c, ("cons_9,"+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		for(const auto& l : sGC.L_t[t]){
			for(const auto& k : sGC.K_l[l]){
				SCIPaddCoefLinear(sGC.scip, c, sGC.varY_lkt[l][k][t],sGC.L[l].energyDemand);
			}
		}
		for(int i=0; i<sGC.d.cardJ; ++i){
			if((sGC.d.rj[i]<=t) && (t<sGC.d.dj[i])) SCIPaddCoefLinear(sGC.scip, c, sGC.varX_it[i][t],-sGC.d.Djk[i][0]);
		}
		SCIPaddCons(sGC.scip, c);
		cons_9.push_back(c);
	}
	sGC.phi_t = vdt;
	sGC.farkas9_t = vdt;
	sGC.cons_9 = cons_9;
	//cout<<"cons_9 ok"<<endl;

	//sGC.nbconstmodele = SCIPgetNConss(sGC.scip);
	//sGC.nbvarmodele = SCIPgetNVars(sGC.scip);
    //sGC.nbcolgenerated = 0;
    
	
	//FILE * filed;
	//filed = fopen("sol_comp", "w");
	//SCIPprintOrigProblem(sGC.scip, filed, "lp", false);



	/*SCIP_CONS * cdelete;
	SCIPcreateConsLinear(sGC.scip, &cdelete, "delete", 0, 0, 0, 17, 17, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[0][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[0][2],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[0][3],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[0][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[1][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[2][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[2][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[2][7],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[2][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[5][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[5][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[6][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[6][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[7][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[8][2],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[8][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete, sGC.varX_it[9][6],1);
	SCIPaddCons(sGC.scip, cdelete);
	cout << "contrainte avec les x=1 ajoutée ! " << endl;*/

	/*SCIP_CONS * cdelete2;
	SCIPcreateConsLinear(sGC.scip, &cdelete2, "delete2", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	//SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[0][8],1);
	//SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[0][9],1);
	//SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[0][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][0],1);
	/*SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][2],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][3],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][7],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[1][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[2][0],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[2][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[2][2],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[2][3],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[2][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][0],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][7],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[3][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][0],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][7],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[4][10],1);* /
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[5][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[5][7],1);
	/*SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[5][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[5][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[5][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][0],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][2],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][3],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[6][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][3],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][6],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][7],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[7][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][0],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][1],1);* /
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][3],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][7],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][8],1);
	/*SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[8][10],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][0],1);* /
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][1],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][4],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][5],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][7],1);
	/*SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][8],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][9],1);
	SCIPaddCoefLinear(sGC.scip, cdelete2, sGC.varX_it[9][10],1);* /
	SCIPaddCons(sGC.scip,cdelete2);
	cout << "contrainte avec les x=0 ajoutée"<<endl;*/



	
    return SCIP_OKAY;
	
}