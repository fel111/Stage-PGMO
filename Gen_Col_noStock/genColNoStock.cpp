#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
#include <limits>
//#include <math.h>
//#include <algorithm>
#include <chrono>
//#include <deque>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "../struct.h"
#include "struct_gencol.h"
#include "../Lecteur_Fichiers/lecteur_taches.h"
#include "../Lecteur_Fichiers/lecteur_param.h"
#include "pricer.h"
#include "../Modele_compact/modele_entier_cplex.h"
#include "../Ordo_compact/ordo_cplex.h"
//#include "genColNoStock.h"

using namespace std;
float infini = numeric_limits<float>::max();

SCIP_RETCODE Load_Original_Model(structGenCol & sGC)
{
	SCIP_PROBDATA* probdata = NULL;
	
	/* create environment */
    SCIP_CALL( SCIPcreate(&sGC.scip));
    assert(sGC.scip != NULL);
    
    /* Load plugin */
    SCIP_CALL( SCIPincludeDefaultPlugins(sGC.scip) );
    if(sGC.p->aff_log_ordo == 0) SCIPsetMessagehdlr(sGC.scip,NULL);
	SCIP_CALL( SCIPchgRealParam(sGC.scip,SCIPgetParam(sGC.scip,"limits/time"),sGC.p->time_limit_ordo) );
	SCIP_CALL( SCIPsetIntParam(sGC.scip, "presolving/maxrounds", 0));
	//SCIP_CALL( SCIPchgRealParam(sGC.scip,SCIPgetParam(sGC.scip,"numerics/epsilon"),0.0001) );
    SCIP_CALL( SCIPchgIntParam(sGC.scip,SCIPgetParam(sGC.scip,"lp/threads"),1) );
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
	for(int i=0; i<sGC.d->cardJ; ++i){
		vector<SCIP_VAR *> v (sGC.d->cardT,NULL);
		for(int t=sGC.d->rj[i]; t<sGC.d->dj[i]; ++t){
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
		vector<SCIP_VAR *> v (sGC.d->cardT, NULL);
		y_lt.push_back(v);
		for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
			SCIP_VAR * varr;
			SCIPcreateVarBasic(sGC.scip, &varr, ("y_lt"+to_string(l)+','+to_string(t)).c_str(), 0.0, SCIPinfinity(sGC.scip), p_t(t,sGC.L[l].energyDemand,sGC), SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(sGC.scip,varr);
			y_lt[l][t] = varr;
		}
	}
	sGC.varY_lt = y_lt;

    

	// Add constraints
    
	
	// x_it - sum_l a_il y_lt = 0
	vector<vector<SCIP_CONS *> > cons_1;
	for(int i=0; i<sGC.d->cardJ; ++i){
		vector<SCIP_CONS *> c;
		cons_1.push_back(c);
		for(int t=sGC.d->rj[i]; t<sGC.d->dj[i]; ++t){
			SCIP_CONS * ct;
			SCIPcreateConsLinear(sGC.scip, &ct, ("cons_1,"+to_string(i)+","+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(sGC.scip, ct, sGC.varX_it[i][t],1);		
			for(const auto& l : sGC.L_t[t]){
				SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lt[l][t],-sGC.a_il[i][l]);
			}
			SCIP_CALL(SCIPaddCons(sGC.scip, ct));
			cons_1[i].push_back(ct);
		}
		vector<double> vd (sGC.d->dj[i]-sGC.d->rj[i],0.0);
		sGC.w_it.push_back(vd);
	}
	sGC.cons_1 = cons_1;

	// sum_t sum_l a_il y_lt >= pi
	vector<SCIP_CONS *> cons_2;
	for(int i=0; i<sGC.d->cardJ; ++i){
		SCIP_CONS * c;
		SCIPcreateConsLinear(sGC.scip, &c, ("cons_2,"+to_string(i)).c_str(), 0, 0, 0, sGC.d->pj[i], SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
		
		for(int l=0; l<sGC.L.size(); ++l){
			for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
				//cout << "a_"<<i<<l<<" : "<<sGC.a_il[i][l]<<endl;
				SCIPaddCoefLinear(sGC.scip, c, sGC.varY_lt[l][t],sGC.a_il[i][l]);
			}
		}
		
		SCIPaddCons(sGC.scip, c);
		cons_2.push_back(c);
	}
	vector<double> vdj (sGC.d->cardJ,0.0);
	sGC.u_i = vdj;
	sGC.cons_2 = cons_2;
	//cout<<"cons_2 ok"<<endl;

	// sum_l -ylt >= -1
	vector<SCIP_CONS *> cons_3;
	for(int t=0; t<sGC.d->cardT; ++t){
		SCIP_CONS * ct;
		SCIPcreateConsLinear(sGC.scip, &ct, ("cons_3,"+to_string(t)).c_str(), 0, 0, 0, -1, SCIPinfinity(sGC.scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  	
		for(const auto& l : sGC.L_t[t]){
			SCIPaddCoefLinear(sGC.scip, ct, sGC.varY_lt[l][t],-1);		
		}
		SCIPaddCons(sGC.scip, ct);
		cons_3.push_back(ct);
	}
	sGC.cons_3 = cons_3;
	vector<double> vdt (sGC.d->cardT,0.0);
	sGC.v_t = vdt; 
	//cout<<"cons_3 ok"<<endl;

	
    return SCIP_OKAY;
	
}

float firstSol(structGenCol & sGC){
    
    //int cptId = 0;
    float upperbound = 0.0;
    
    for(int t=sGC.d->releaseDateMin; t<sGC.d->cardT; ++t){
        vector<int> taskList;
        feasibleSet l;
        for(int i=0; i<sGC.d->cardJ; ++i){
            if((sGC.d->rj[i]<=t) && (t<sGC.d->rj[i]+sGC.d->pj[i])){
                taskList.push_back(i);
                if(taskList.size() == 1){
                    l.energyDemand = sGC.d->Dj[i];
                    l.deadLine = sGC.d->dj[i];
                    l.releaseTime = sGC.d->rj[i];
                }
                else{
                    l.energyDemand += sGC.d->Dj[i];
                    if(l.deadLine > sGC.d->dj[i]) l.deadLine = sGC.d->dj[i];
                    if(l.releaseTime < sGC.d->rj[i]) l.releaseTime = sGC.d->rj[i];
                }
            }
        }
        if(taskList.size() > 0){
            //l.id = cptId;
            //l.timeGen = t; 
            l.tasksList = taskList;
            if(checkSet(l,sGC)==-1){
                //l.cost = p_t(t,l.energyDemand,sGC);
                //cout << "l, cout : " << l.id <<"  "<<l.cost << endl;
                l.id = sGC.L.size();
                sGC.L.push_back(l);
                sGC.cardL = sGC.L.size();
                //addSetK_l(l,sGC);
                addA_il(l,sGC);
                addL_t(l,sGC);
                //cptId++;
            }//else cout << "set deja present " << l.id << endl;
            upperbound += p_t(t,l.energyDemand,sGC);
        }
	}
    return upperbound;
}


/*int initData(structGenCol & sGC, string nb){
    data d;
    param p;

    d.cardR = 2;
	d.s0 = 0;
	d.cardM = 1;
    string param = "param_boucle_ordo_ls.txt";
    //string instance = "inst_100";
	if(lecteur_param("Param/"+param,p,d) == 0) return 0;
	//d.Q = 20;
	//p.qmax = 20;
	//p.qmin = 0;
	//p.qinit = 0;
	
	lecteur_taches_EnergSchedInst("Donnees/dataSchedulingInstances_newname/inst_"+nb,d);

    for(int i=0;i<d.cardT;++i){
        
		/*vector<float> temp_bpt = {0.0,4.0,4.0,8.0,8.0,100.0,100.0,infini};
		vector<float> temp_valbpt = {0.0,4.0,4.0,6.0,6.0,190.0};* /
		vector<float> temp_pente = {1.0,0.5,2.0,1.5};

		vector<float> temp_bpt = {0.0,4.0,8.0,100.0,infini};
		vector<float> temp_valbpt = {0.0,4.0,6.0,190.0};
		//vector<float> temp_bpt = {0.0,5.0,infini};
		//vector<float> temp_valbpt = {0.0,10.0};
		//vector<float> temp_pente = {2.0,1.5};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(4);
	}
	//vector<vector<float> > ress;
	//for(int i=0;i<d.cardM;++i){
		//vector<float> r = {(d.cardJ*4.0/d.cardM),(d.cardJ*4000.0/d.cardM)}; // s'assure qu'il y ai assez de ram/cpu
		//vector<float> r = {d.cardJ, d.cardJ};
		//cout << r[0] << " " << r[1] << endl;
		//ress.push_back(r);
	//}
	//d.Ckr = ress;
	//vector<float> Dk (d.cardM, 0.0); // ne converge pas avec != 0
	//d.Dk = Dk;

    sGC.d = d;
    sGC.p = p;

    

    //initialisatin P_0
    //addP_0(sGC);

}*/


float genColNoStock(data & d,param & p,float &tps,vector<float> & demande, float &gap, string &statut){
    structGenCol sGC(d,p);
    //initialisation L_t
    for(int t=0; t<sGC.d->cardT; ++t){
        vector<int> lt_temp;
        sGC.L_t.push_back(lt_temp);
    }
    //initialisation a_il
    for(int j=0; j<sGC.d->cardJ; ++j){
        vector<int> ailtemp;
        sGC.a_il.push_back(ailtemp);
    }

    // affichage fonction pwl
    /*for(int t=0; t<d.cardT; ++t){
        cout << "bpt, valbpt " << t << endl;
        for(int i=0; i<d.nb_bp[t]; ++i){
            cout << d.bpt[t][i] << " ";
            cout << d.valbpt[t][i] << " ";
        }
        cout << endl;
    }*/

	int tmp;
    SCIP_RETCODE status;
    SCIP_STATUS solstatus;
    //initData(sGC);
    
    //cout <<"initData OK ! "<<endl;
	
	float bsup = firstSol(sGC);
    //cout << "bsup firstSol : " << bsup << endl;
	//cout <<"firstSol OK ! "<<endl;
	
	//affK_l(sGC);
	//affL_t(sGC);
	//affAllSet(sGC);
	


    // Load and build model
    SCIP_CALL (Load_Original_Model(sGC));
    
    // Set dedicated parameters
    //SCIP_CALL (status = SetColGenParameters_Scip((*pU).scip, (*pU).Params));
    
    
    // activate the pricer
    SCIP_CALL( SCIPactivatePricer(sGC.scip, SCIPfindPricer(sGC.scip, "generatetestnewobjects") ));
    
    // ajout borne sup solution initiale
    SCIPsetObjlimit(sGC.scip,bsup);

    // SOLVE
    //cout << "Start solve..."<<endl;
    auto start_time = chrono::steady_clock::now();
    SCIP_CALL( status = SCIPsolve(sGC.scip) );
    auto end_time = chrono::steady_clock::now();
    tps = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;


    sGC.sol = SCIPgetBestSol(sGC.scip);

    //cout << "verif : "<<verifSol(sGC) << endl; 

    
    //cout << "obj : "<< SCIPgetSolOrigObj(sGC.scip,sGC.sol) << endl;
    /*for(int j=0; j<sGC.d->cardJ; ++j){
        for(int t=sGC.d->rj[j]; t<sGC.d->dj[j]; ++t){
            if(SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][t]) >= 0.9) cout <<"tache "<<j<<" réalisée au temps "<<t<<endl;
        }
    }*/
    
    
    // export modele dernier PMR
    //SCIPwriteTransProblem(sGC.scip,"lastPMRNoStock.lp","lp",0);
    /*FILE * filed;
	filed = fopen("afterCol", "a+");
	SCIPprintTransProblem(sGC.scip, NULL, "lp", 0);
    fclose(filed);*/


    //cout << "------sol SCIP------"<<endl;
    // aff solution
    //SCIPprintSol(sGC.scip, sGC.sol, NULL, 0);

    // Get the current scip solving time in seconds.
    //cout << "temps scip : "<<SCIPgetSolvingTime(sGC.scip) <<endl;
    
    
    //affAllSet(sGC);

    // export du dernier pb maitre
    /*FILE * filed;
    filed = fopen("lastPM","w");
    SCIPprintTransProblem(sGC.scip, filed, NULL, 0);
    fclose(filed);*/

    //cout << "---------" << endl;

    /*for(int l=0; l<sGC.L.size(); ++l){
        for(const auto& k : sGC.K_l[l]){
		    for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
	            if(sGC.varY_lkt[l][k][t] != NULL){
                    //if(SCIPisGT(sGC.scip,SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY_lkt[l][k][t]),0) ) cout << "y"<<l<<","<<k<<","<<t<<" = "<<SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY_lkt[l][k][t])<<endl;
                    cout << "y"<<l<<","<<k<<","<<t<<" = "<<SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY_lkt[l][k][t])<<endl;
                }
            }
		}
	}
    cout << "---------" << endl;
    for(int k=0; k<sGC.d->nb_bp[0]; ++k){
		for(int t=0; t<sGC.d->cardT; ++t){
            if(SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY0_kt[k][t])>0) cout << "y0,"<<k<<","<<t<<" = "<<SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY0_kt[k][t])<<endl;
        }
    }*/

    //SCIPprintStatistics(sGC.scip,NULL);

    /*cout << "dual bound : " << SCIPgetDualbound(sGC.scip);
    cout << "primal bound : " << SCIPgetPrimalbound(sGC.scip);
    cout << "scip time : " << SCIPgetSolvingTime(sGC.scip);
    cout << "npricevar : " << SCIPgetNPricevarsFound(sGC.scip);



    cout << "nb node : " <<SCIPgetNTotalNodes(sGC.scip)<<endl;
    cout << "nb col gen : " << sGC.nbcolgenerated << endl;
    cout << "cptPricer : "<<sGC.cptPricer<<endl;
    cout << "temps gencol : "<<tpsgencol<<endl;
    cout << "temps scip : "<<SCIPgetSolvingTime(sGC.scip) <<endl;
    cout << "temps pricer : "<<sGC.tpsPricer<<endl;
    cout << "temps SP : "<<sGC.tpsSP<<endl;*/


    // release variables
    /*int fois = 0;
    for(int i=0; i<sGC.d->cardJ; ++i){
        for(int t=0; t<sGC.d->cardT; ++t){
            if((sGC.d->rj[i]<=t)&&(t<sGC.d->dj[i])){
                SCIPreleaseVar(sGC.scip,&sGC.varX_it[i][t]);
                SCIPreleaseCons(sGC.scip,&sGC.cons_1[i][t-sGC.d->rj[i]]);
            }
            if(fois==0) SCIPreleaseCons(sGC.scip,&sGC.cons_3[t]);
        }
        SCIPreleaseCons(sGC.scip,&sGC.cons_2[i]);
        ++fois;
    }
    for(int l=0; l<sGC.L.size(); ++l){
		for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
            if(sGC.varY_lt[l][t] != NULL) SCIPreleaseVar(sGC.scip,&sGC.varY_lt[l][t]);
        }
    }*/
    float obj = -1.0;
    if (status == SCIP_OKAY){
        if(SCIPgetStatus(sGC.scip) == SCIP_STATUS_OPTIMAL) statut = "optimal";
        if(SCIPgetStatus(sGC.scip) == SCIP_STATUS_TIMELIMIT) statut = "feasible";
        gap = SCIPgetGap(sGC.scip);
        obj = SCIPgetSolOrigObj(sGC.scip,sGC.sol);
        vector<float> prod;
        for(int t=0; t<sGC.d->cardT; ++t){
            float dem=0.0;
            for(int j=0; j<sGC.d->cardJ; ++j){
                if(sGC.varX_it[j][t] != NULL){
                    if(SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][t]) >= 0.9) dem += sGC.d->Dj[j];
                }
            }
            prod.push_back(dem);
            //cout << "gencol : demande " << prod[t] << endl;
        }
        demande = prod;
    }
    SCIP_CALL ( SCIPfree(&sGC.scip) );
    return obj;
}




/*int main(int argc, char* argv[]){
    string nb = argv[1];
    float tpsGenCol;
    float solGenCol;
    vector<float> temp;
    float gapGenCol;
    string statutGenCol;
    structGenCol sGC;
    initData(sGC,nb);
    solGenCol = genColNoStock(sGC.d,sGC.p,tpsGenCol,temp,gapGenCol,statutGenCol);
    cout << solGenCol << " " << tpsGenCol << " " << statutGenCol << " " << gapGenCol << " // GENCOL : solution, temps, statut, gap" << endl;

    float tpsComp, solComp, gapComp;
    string statutComp;
    solComp = ordo_cplex(sGC.d,sGC.p,tpsComp,statutComp,gapComp);
    cout << solComp << " " << tpsComp << " " << statutComp << " " << gapComp << " // COMPACT : solution, temps, statut, gap" << endl;


    return 0;
}*/