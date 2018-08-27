#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
#include <chrono>
//#include <limits>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "../struct.h"
//#include "lotsizcontcom.h"
#include <ilcplex/ilocplex.h>
#include "ordo_cplex2.h"

ILOSTLBEGIN


float ordo_cplex(data const& d,param const& p, float &tps, vector<float> &demande,string &status, float &gap){
	IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    if (p.aff_log_compact==0) cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads,p.nb_threads_cplex);
    cplex.setParam(IloCplex::TiLim,p.time_limit_ordo);

	// affichage fonction pwl
    /*for(int t=0; t<d.cardT; ++t){
        cout << "bpt, valbpt " << t << endl;
        for(int i=0; i<d.nb_bp[t]; ++i){
            cout << d.bpt[t][i] << " ";
            cout << d.valbpt[t][i] << " ";
        }
        cout << endl;
    }*/


	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	//ajout variables xt
	//IloNumVarArray xt (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

/*
	//ajout variables y_jkt binaires
	vector<vector<IloNumVarArray> > y_jkt;
	for(int j=0; j<d.cardJ; ++j){
		vector<IloNumVarArray> v;
		y_jkt.push_back(v);
		for(int k=0; k<d.cardM; ++k){//"yk"+to_string(j)+"_"+to_string(k)+"_"
			IloNumVarArray var (env,d.cardT,0,1,ILOBOOL );
			y_jkt[j].push_back(var);
		}	
	}

	//ajout variables z_kt binaires
	vector<IloNumVarArray> z_kt;
	for(int k=0; k<d.cardM; ++k){// , "zkt"+to_string(k)+"_"
		IloNumVarArray var (env,d.cardT,0,1,ILOBOOL);
		z_kt.push_back(var);
	}
*/

	//ajout variables y_jt binaires
	vector<IloNumVarArray> y_jt;
	for(int j=0; j<d.cardJ; ++j){
		IloNumVarArray var (env,d.cardT,0,1,ILOBOOL );
		y_jt.push_back(var);
	}



	// CONTRAINTES

	//contraintes cout
	for(int t=0; t<d.cardT; ++t){
		IloExpr SumDem(env);
		for(int j=0; j<d.cardJ; ++j){
			SumDem += y_jt[j][t]*d.Dj[j];
		}
		IloNumArray bpt(env, d.nb_bp[t]-1);
        for(int i=0; i<(d.nb_bp[t]-1);++i){
            bpt[i] = d.bpt[t][i+1];
        }

        IloNumArray pente(env, d.nb_bp[t]);
        for(int i=0; i<d.nb_bp[t];++i){
            pente[i] = d.pente[t][i];
        }
		/*if(d.bpt[t][0] > 0){
			IloNumArray pente(env, d.nb_bp[t]+1);
			IloNumArray bpt(env, d.nb_bp[t]);
			pente[0] = 0.0;
			bpt[0] = d.bpt[t][0];
			for(int j=1; j<d.nb_bp[t]+1;++j){
        		pente[j] = d.pente[t][j-1];
				if (j<d.nb_bp[t]) bpt[j] = d.bpt[t][j];
    		}
			model.add(ct[t] == IloPiecewiseLinear(SumDem,bpt,pente,0, 0));
		}
		else{
			IloNumArray pente(env, d.nb_bp[t]);
			IloNumArray bpt(env, d.nb_bp[t]-1);
			for(int j=0; j<d.nb_bp[t];++j){
       			pente[j] = d.pente[t][j];
				if(j<d.nb_bp[t]-1) bpt[j] = d.bpt[t][j+1];
    		}
			model.add(ct[t] == IloPiecewiseLinear(SumDem,bpt,pente,d.bpt[t][0], d.valbpt[t][0]));
		}*/
		model.add(ct[t] == IloPiecewiseLinear(SumDem,bpt,pente,d.bpt[t][0], 0));
	}

    IloExpr obj(env);
    obj = IloSum(ct);
    model.add(IloMinimize(env, obj));
    obj.end();
/*
	// sum y_jkt <= 1
	for(int j=0; j<d.cardJ; ++j){
		for(int t=0; t<d.cardT; ++t){
            IloExpr sumy(env);
            for(int k=0; k<d.cardM; ++k){
                sumy += y_jkt[j][k][t];
            }
            model.add(sumy <= 1);
        }
	}
*/


	// sum_t y_jt >= pj
	for(int j=0; j<d.cardJ; ++j){
        IloExpr sumy2(env);
		for(int t=0; t<d.cardT; ++t){
			sumy2 += y_jt[j][t];
		}		
		model.add(sumy2 >= d.pj[j]);
	}
if(p.taches_avec_fenetre_temps == 1){
	// contrainte respect fenetre temps
	for(int j=0; j<d.cardJ; ++j){
        IloExpr s1(env);
		IloExpr s2(env);
		for(int t=0; t<d.rj[j]; ++t){
			s1 += y_jt[j][t];
		}
		for(int t=d.dj[j]; t<d.cardT; ++t){
			s2 += y_jt[j][t];
		}			
		model.add(s1 + s2 == 0);
	}
}
/*
	// sum_j cjr*yjkt <= Ckr
	for(int k=0; k<d.cardM; ++k){
		for(int r=0; r<d.cardR; ++r){
			for(int t=0; t<d.cardT; ++t){
				IloExpr sumjr(env);
                for(int j=0; j<d.cardJ; ++j){
					sumjr += d.cjr[j][r]*y_jkt[j][k][t];
				}
				model.add(sumjr <= d.Ckr[k][r]);
			}
		}
	}

	// zkt - yjkt >= 0
	for(int k=0; k<d.cardM; ++k){
		for(int j=0; j<d.cardJ; ++j){
			for(int t=0; t<d.cardT; ++t){
                model.add(z_kt[k][t] - y_jkt[j][k][t] >= 0);
			}
		}
	}

	// xt - sum_j sum_k Dj*yjkt - sum_k Dk*zkt = 0
	for(int t=0; t<d.cardT; ++t){
		IloExpr sum1(env);
        for(int j=0; j<d.cardJ; ++j){
			for(int k=0; k<d.cardM; ++k){
				sum1 += d.Dj[j][k]*y_jkt[j][k][t];
			}
		}
        IloExpr sum2(env);
		for(int k=0; k<d.cardM; ++k){
			sum2 += d.Dk[k]*z_kt[k][t];
		}			
		model.add(xt[t]-sum1-sum2 == 0);
	}
*/

	auto start_time = chrono::steady_clock::now();
    cplex.solve();
	auto end_time = chrono::steady_clock::now();
	tps = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	float sol;
	if(cplex.getStatus() == IloAlgorithm::Feasible){
		status = "feasible";
		vector<float> prod;
		//cout << "demande : ";	
		for(int t=0; t<d.cardT; ++t){
			float dem=0.0;
			//cout << SCIPgetSolVal(scip,sol,xt[i]) << ", ";
			//prod.push_back(static_cast<int>(ceil(SCIPgetSolVal(scip,sol,xt[i]))));
			//prod.push_back(roundd(SCIPgetSolVal(scip,sol,xt[i])));
			//prod.push_back(roundd(cplex.getValue(xt[i]),5));
			for(int j=0; j<d.cardJ; ++j){
				dem += cplex.getValue(y_jt[j][t])*d.Dj[j];
			}
			prod.push_back(dem);
			//cout << "cplex : demande " << prod[t] << endl;
		}
		demande = prod;
		sol = cplex.getObjValue();
		gap = cplex.getMIPRelativeGap();
	}
	else if(cplex.getStatus() == IloAlgorithm::Optimal){
		status = "optimal";
		vector<float> prod;
		/*//cout << "demande : ";	
		//for(int i=0; i<d.cardT; ++i){
			//cout << SCIPgetSolVal(scip,sol,xt[i]) << ", ";
			//prod.push_back(static_cast<int>(ceil(SCIPgetSolVal(scip,sol,xt[i]))));
			//prod.push_back(roundd(SCIPgetSolVal(scip,sol,xt[i])));
			//prod.push_back(roundd(cplex.getValue(xt[i]),5));
			//  prod.push_back(cplex.getValue(xt[i]));
			//cout << prod[i] << " ";
			//cout << "------------ CPLEX SOL ----------"<<endl;
			for(int t=0; t<d.cardT; ++t){
				//if(cplex.getValue(xt[t])> 0) cout << "xt" << t << ": " <<cplex.getValue(xt[t]) << endl;
				//if(cplex.getValue(ct[t]) > 0) cout << "ct" << t << ": " <<cplex.getValue(ct[t]) << endl;
				//if(cplex.getValue(st[t]) > 0 ) cout << "st" << t << ": " <<cplex.getValue(st[t]) << endl;
				for(int k=0; k<d.cardM; ++k){
					//if(cplex.getValue(z_kt[k][t]) > 0.5) cout << "zkt" << k << t << ": " <<cplex.getValue(z_kt[k][t]) << endl;
					for(int j=0; j<d.cardJ; ++j){
						if(cplex.getValue(y_jkt[j][k][t]) > 0.9) cout << "yjkt" << j << k<<t << ": " <<cplex.getValue(y_jkt[j][k][t]) << endl;
					}
				}
			}
		//} */
		for(int t=0; t<d.cardT; ++t){
			float dem=0.0;
			//cout << SCIPgetSolVal(scip,sol,xt[i]) << ", ";
			//prod.push_back(static_cast<int>(ceil(SCIPgetSolVal(scip,sol,xt[i]))));
			//prod.push_back(roundd(SCIPgetSolVal(scip,sol,xt[i])));
			//prod.push_back(roundd(cplex.getValue(xt[i]),5));
			for(int j=0; j<d.cardJ; ++j){
				dem += cplex.getValue(y_jt[j][t])*d.Dj[j];
			}
			prod.push_back(dem);
			//cout << "cplex : demande " << prod[t] << endl;
			//cout << prod[i] << " ";
		}
		demande = prod;
		sol = cplex.getObjValue();
		gap = cplex.getMIPRelativeGap();
		//cout << "nb noeuds cplex : " << cplex.getNnodes() << endl;
	}
	else{
		status = "Unknown";
		sol = -1.0;
	}
	env.end();
	
	return sol;


}