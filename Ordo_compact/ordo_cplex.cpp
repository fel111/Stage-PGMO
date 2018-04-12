#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
//#include <limits>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "data_struct.h"
//#include "lotsizcontcom.h"
#include <ilcplex/ilocplex.h>
#include "ordo_cplex.h"

ILOSTLBEGIN


float ordo_cplex(data &d){
	IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    //cplex.setParam(IloCplex::Threads, 1);
    cplex.setParam(IloCplex::TiLim,300);


	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	//ajout variables xt
	IloNumVarArray xt (env,d.cardT,0,IloInfinity,ILOFLOAT);


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
	for(int i=0; i<d.cardT; ++i){
		if(d.bpt[i][0] > 0){
			IloNumArray pente(env, d.nb_bp[i]+1);
			IloNumArray bpt(env, d.nb_bp[i]);
			pente[0] = 0.0;
			bpt[0] = d.bpt[i][0];
			for(int j=1; j<d.nb_bp[i]+1;++j){
        		pente[j] = d.pente[i][j-1];
				if (j<d.nb_bp[i]) bpt[j] = d.bpt[i][j];
    		}
			model.add(ct[i] == IloPiecewiseLinear(xt[i],bpt,pente,0, 0));
		}
		else{
			IloNumArray pente(env, d.nb_bp[i]);
			IloNumArray bpt(env, d.nb_bp[i]-1);
			for(int j=0; j<d.nb_bp[i];++j){
       			pente[j] = d.pente[i][j];
				if(j<d.nb_bp[i]-1) bpt[j] = d.bpt[i][j+1];
    		}
			model.add(ct[i] == IloPiecewiseLinear(xt[i],bpt,pente,d.bpt[i][0], d.valbpt[i][0]));
		}
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

	// xt - sum_j sum_k Djk*yjkt - sum_k Dk*zkt = 0
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
		model.add(xt[t]-sum1-sum2 == 0);
	}

    cplex.solve();

	vector<float> prod;
	cout << "demande : ";	
	for(int i=0; i<d.cardT; ++i){
		//cout << SCIPgetSolVal(scip,sol,xt[i]) << ", ";
		//prod.push_back(static_cast<int>(ceil(SCIPgetSolVal(scip,sol,xt[i]))));
		//prod.push_back(roundd(SCIPgetSolVal(scip,sol,xt[i])));
		prod.push_back(roundd(cplex.getValue(xt[i]),5));
        cout << prod[i] << " ";
	}
	cout << endl;

	d.dt = prod;

    float sol = cplex.getObjValue();
    //env.out() << "Cost LotsizEnt = " << cplex.getObjValue() << endl;

    env.end();
    return sol;

	

}
