#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <chrono>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "data_struct.h"
//#include "lotsizcontcom.h"
#include <ilcplex/ilocplex.h>
#include "lotsizentCPX.h"

//#define SCIP_DEBUG
ILOSTLBEGIN

float lotsizentCPX(data const& d, param const& p, vector<float> &varia, float &tps, float &borneinf, string &status){
    vector<int> dt = dtToInt(d,p.choix_dt_ls);
    int batMax;
    if(p.choix_dt_ls == 2) batMax = d.Q*100;
    else batMax = d.Q;

    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads,p.nb_threads_cplex);
    cplex.setParam(IloCplex::TiLim,p.time_limit_lotsizingent);
   
	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);
	
	//ajout variables st de stock
    IloNumVarArray st (env,d.cardT,0,batMax,ILOINT);

	//ajout variables xt de production
	IloNumVarArray xt (env,d.cardT,0,IloInfinity,ILOINT);

	
	//CONTRAINTES
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
	
	//contraintes bilan energetique
	model.add(xt[0] + d.s0 - st[0] == dt[0]);
	for(int i=1; i<d.cardT; ++i){
         model.add(xt[i] + st[i-1] - st[i] == dt[i]);
	}


	//contrainte stmax  == 0
	model.add(st[d.cardT-1] == d.s0);


    auto start_time = chrono::steady_clock::now();
    cplex.solve();
	auto end_time = chrono::steady_clock::now();
	tps = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;


    float sol;
	if(cplex.getStatus() == IloAlgorithm::Feasible){
		status = "Feasible";
		borneinf = cplex.getBestObjValue();
		sol = cplex.getObjValue();
		vector<float> variat;
        if(p.choix_dt_ls != 2){
            variat.push_back(roundd(cplex.getValue(st[0]),5) - d.s0);
            for(int t=1; t<d.cardT; ++t){
                variat.push_back(roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5));
            }
            varia = variat;
        }
        else{
            variat.push_back((roundd(cplex.getValue(st[0]),5) - d.s0)/100);
            for(int t=1; t<d.cardT; ++t){
                variat.push_back(roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5)/100);
            }
            varia = variat;
        }
	}
    else if(cplex.getStatus() == IloAlgorithm::Optimal){
		status = "Optimal";
		borneinf = cplex.getBestObjValue();
		sol = cplex.getObjValue();
		vector<float> variat;
        if(p.choix_dt_ls != 2){
            variat.push_back(roundd(cplex.getValue(st[0]),5) - d.s0);
            for(int t=1; t<d.cardT; ++t){
                variat.push_back(roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5));
            }
            varia = variat;
        }
        else{
            variat.push_back((roundd(cplex.getValue(st[0]),5) - d.s0)/100);
            for(int t=1; t<d.cardT; ++t){
                variat.push_back(roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5)/100);
            }
            varia = variat;
        }
	}
	else{
		status = "Unknown";
		borneinf = cplex.getBestObjValue();
		sol = -1.0;
	}

	
    env.end();
	return sol;

}












/*
float lotsizentCPX(data const& d, int choix,vector<float> & var, param const& p){
    vector<int> dt = dtToInt(d,choix);

    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads,p.nb_threads_cplex);
    cplex.setParam(IloCplex::TiLim,p.time_limit_lotsizingent);
   
	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);
	
	//ajout variables st de stock
    IloNumVarArray st (env,d.cardT,0,d.Q,ILOINT);

	//ajout variables xt de production
	IloNumVarArray xt (env,d.cardT,0,IloInfinity,ILOINT);

	
	//CONTRAINTES
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
	
	//contraintes bilan energetique
	model.add(xt[0] +d.s0 - st[0] == dt[0]);
	for(int i=1; i<d.cardT; ++i){
         model.add(xt[i] + st[i-1] - st[i] == dt[i]);
	}


	//contrainte stmax  == 0
	model.add(st[d.cardT-1] == d.s0);

		
    //IloCplex cplex(env);
    //cplex.extract(model);
    //cplex.exportModel("cplex_lotsizent.lp");

    //auto start_time = chrono::steady_clock::now();
    cplex.solve();
    //auto end_time = chrono::steady_clock::now();
    //cout << "temps com : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;

    //cout << "nb noeuds explorés: "<< cplex.getNnodes()<<endl;
    //borninf = cplex.getBestObjValue();
    //status = cplex.getStatus();

    vector<float> variat (d.cardT,0.0);
	
	variat[0] = roundd(cplex.getValue(st[0]),5) - d.s0;
	for(int t=1; t<d.cardT; ++t){
		variat[t] = roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5);
	}
	var = variat;

    float sol = cplex.getObjValue();
    //env.out() << "Cost LotsizEnt = " << cplex.getObjValue() << endl;

    env.end();
    return sol;

}*/