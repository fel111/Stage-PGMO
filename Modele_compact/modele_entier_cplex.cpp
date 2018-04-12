#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "data_struct.h"
//#include "compact.h"
#include <ilcplex/ilocplex.h>
//#define SCIP_DEBUG
using namespace std;
//int infini_int = numeric_limits<int>::max();


void modele_entier_cplex(data d){

    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads, 1);
    cplex.setParam(IloCplex::TiLim,300);


	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	//ajout variables xt
	IloNumVarArray xt (env,d.cardT,0,IloInfinity,ILOFLOAT);

    //ajout variables st de stock
    IloNumVarArray st (env,d.cardT,0.0,d.Q,ILOFLOAT);


	//ajout variables y_jkt binaires
	vector<vector<IloNumVarArray> > y_jkt;
	for(int j=0; j<d.cardJ; ++j){
		vector<IloNumVarArray> v;
		y_jkt.push_back(v);
		for(int k=0; k<d.cardM; ++k){
			IloNumVarArray var (env,d.cardT,0,1,ILOBOOL);
			y_jkt[j].push_back(var);
		}	
	}

	//ajout variables z_kt binaires
	vector<IloNumVarArray> z_kt;
	for(int k=0; k<d.cardM; ++k){
		IloNumVarArray var (env,d.cardT,0,1,ILOBOOL);
		z_kt.push_back(var);
	}

	// CONTRAINTES

	//contraintes sur ct
	IloNumArray bpt(env, d.nb_bp[0]-1);
    for(int i=0; i<(d.nb_bp[0]-1);++i){
        bpt[i] = d.bpt[0][i+1];
    }

    IloNumArray pente(env, d.nb_bp[0]);
    for(int i=0; i<d.nb_bp[0];++i){
        pente[i] = d.pente[0][i];
    }

	for(int i=0; i<d.cardT; ++i){
        model.add(ct[i] == IloPiecewiseLinear(xt[i],
                                       bpt,
                                       pente,
                                       0.0, 0.0));
	}

    IloExpr obj(env);
    obj = IloSum(ct);
    model.add(IloMinimize(env, obj));
    obj.end();

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

	// sum_k sum_t y_jkt = pj
	for(int j=0; j<d.cardJ; ++j){
        IloExpr sumy2(env);
        for(int k=0; k<d.cardM; ++k){
			for(int t=0; t<d.cardT; ++t){
				sumy2 += y_jkt[j][k][t];
			}
		}			
		model.add(sumy2 == d.pj[j]);
	}

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

	// xt - +st -st-1 - sum_j sum_k Djk*yjkt - sum_k Dk*zkt = 0
	
    
    for(int t=0; t<d.cardT; ++t){
		IloExpr sum1(env);
        for(int j=0; j<d.cardJ; ++j){
			for(int k=0; k<d.cardM; ++k){
				sum1 += d.Djk[j][k]*y_jkt[j][k][t];
			}
		}
        IloExpr sum2(env);
		for(int k=0; k<d.cardM; ++k){
			sum2 += d.Dk[k]*z_kt[k][t];
		}
        if(t!=0) model.add(xt[t]+st[t]-st[t-1]-sum1-sum2 == 0);
        else model.add(xt[t]+st[t]-d.s0-sum1-sum2 == 0);
	}

    model.add(st[d.cardT-1]-d.s0 >= 0);

    cplex.solve();

    cout << "obj modele entier cplex : "<<cplex.getObjValue()<<endl;
    //env.out() << "Cost LotsizEnt = " << cplex.getObjValue() << endl;

    env.end();
    
}