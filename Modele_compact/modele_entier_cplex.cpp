#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <chrono>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "../struct.h"
//#include "compact.h"
#include <ilcplex/ilocplex.h>
//#define SCIP_DEBUG
using namespace std;
//int infini_int = numeric_limits<int>::max();


float modele_entier_cplex(data const& d, param const& p, float &tps, float &borneinf, string &status){

    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    if (p.aff_log_compact_cplex==0) cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads,p.nb_threads_cplex);
    cplex.setParam(IloCplex::TiLim,p.time_limit_compact);
	//cplex.setParam(IloCplex::NodeLim,0);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve,0);


	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	//ajout variables xt
	IloNumVarArray xt (env,d.cardT,0,IloInfinity,ILOFLOAT);

    //ajout variables st de stock
    IloNumVarArray st (env,d.cardT+1,0.0,d.Q,ILOFLOAT);


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
if(p.taches_avec_fenetre_temps == 1){
	// contrainte respect fenetre temps
	for(int j=0; j<d.cardJ; ++j){
        IloExpr s1(env);
		IloExpr s2(env);
        for(int k=0; k<d.cardM; ++k){
			for(int t=0; t<d.rj[j]; ++t){
				s1 += y_jkt[j][k][t];
			}
			for(int t=d.dj[j]; t<d.cardT; ++t){
				s2 += y_jkt[j][k][t];
			}
		}			
		model.add(s1 + s2 == 0);
	}
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

	// xt -st +st-1 - sum_j sum_k Djk*yjkt - sum_k Dk*zkt = 0
	/*
    
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
        if(t!=0) model.add(xt[t]-st[t]+st[t-1]-sum1-sum2 == 0);
        else model.add(xt[t]-st[t]+d.s0-sum1-sum2 == 0);
	}*/

	// st+1 - st - xt + conso
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
        //if(t!=0) model.add(xt[t]-st[t]+st[t-1]-sum1-sum2 == 0);
        //else model.add(xt[t]-st[t]+d.s0-sum1-sum2 == 0);
		model.add(st[t+1]-st[t]-xt[t]+sum1+sum2 == 0);
	}

	model.add(st[0] == d.s0);
    model.add(st[d.cardT]-st[0] >= 0);
	//


	//cout << "s0 = "<<d.s0<<endl;
	//cout << "qmax = "<<d.Q<<endl;
	auto start_time = chrono::steady_clock::now();
    cplex.solve();
	auto end_time = chrono::steady_clock::now();
	tps = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	
	float sol;
	if(cplex.getStatus() == IloAlgorithm::Feasible){
		status = "Feasible";
		borneinf = cplex.getBestObjValue();
		sol = cplex.getObjValue();
	}
    else if(cplex.getStatus() == IloAlgorithm::Optimal){
		status = "Optimal";
		borneinf = cplex.getBestObjValue();
		sol = cplex.getObjValue();
		//cout << "------------ CPLEX SOL ----------"<<endl;
		/*for(int t=0; t<d.cardT; ++t){
			if(cplex.getValue(xt[t])> 0) cout << "xt" << t << ": " <<cplex.getValue(xt[t]) << endl;
			if(cplex.getValue(ct[t]) > 0) cout << "ct" << t << ": " <<cplex.getValue(ct[t]) << endl;
			if(cplex.getValue(st[t]) > 0 ) cout << "st" << t << ": " <<cplex.getValue(st[t]) << endl;
			for(int k=0; k<d.cardM; ++k){
				if(cplex.getValue(z_kt[k][t]) > 0.5) cout << "zkt" << k << t << ": " <<cplex.getValue(z_kt[k][t]) << endl;
				for(int j=0; j<d.cardJ; ++j){
					if(cplex.getValue(y_jkt[j][k][t]) > 0.5) cout << "yjkt" << j << k<<t << ": " <<cplex.getValue(y_jkt[j][k][t]) << endl;
				}
			}
		}*/
		/*for(int t=6; t<8; ++t){
			cout << "xt" << t << ": " <<cplex.getValue(xt[t]) << endl;
			cout << "ct" << t << ": " <<cplex.getValue(ct[t]) << endl;
			cout << "st" << t << ": " <<cplex.getValue(st[t]) << endl;
			for(int k=0; k<d.cardM; ++k){
				cout << "zkt" << k << t << ": " <<cplex.getValue(z_kt[k][t]) << endl;
				for(int j=0; j<d.cardJ; ++j){
					cout << "yjkt" << j << k<<t << ": " <<cplex.getValue(y_jkt[j][k][t]) << endl;
				}
			}
		}*/


	}
	else{
		status = "Unknown";
		borneinf = cplex.getBestObjValue();
		sol = -1.0;
	}

	//cplex.exportModel("cplexmodelecompact.lp");
	//cout << "borne inf cplex : " << cplex.getObjValue() <<endl;
    cout << "nb noeuds cplex (compact) : " << cplex.getNnodes() << endl;
	env.end();

    return sol;
}












float relaxation_modele_entier_cplex(data const& d,vector<float>& relax, param const& p, float &tps, float &borneinf, string &status){

    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    cplex.setOut(env.getNullStream());
    cplex.setParam(IloCplex::Threads,p.nb_threads_cplex);
    cplex.setParam(IloCplex::TiLim,p.time_limit_relax);
	//cplex.setParam(IloCplex::NodeLim,0);
	//cplex.setParam(IloCplex::Param::Preprocessing::Presolve,0);
	//cplex.setParam(IloCplex::IntSolLim,1);


	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	//ajout variables xt
	IloNumVarArray xt (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	//variable xtij pour les fcts par morceaux
	vector<IloNumVarArray> xij;
	for(int i=0; i<d.cardT; ++i){
		IloNumVarArray var (env,d.nb_bp[i],0.0,IloInfinity,ILOFLOAT);
		xij.push_back(var);
	}

    //ajout variables st de stock
    IloNumVarArray st (env,d.cardT,0.0,d.Q,ILOFLOAT);

	//ajout variables binaires pwd
	vector<IloNumVarArray>  pwd;
	for(int i=0; i<d.cardT; ++i){
		IloNumVarArray var (env,d.nb_bp[i],0.0,1.0,ILOBOOL);
		pwd.push_back(var);
	}

	//ajout variables y_jkt binaires
	vector<vector<IloNumVarArray> > y_jkt;
	for(int j=0; j<d.cardJ; ++j){
		vector<IloNumVarArray> v;
		y_jkt.push_back(v);
		for(int k=0; k<d.cardM; ++k){
			IloNumVarArray var (env,d.cardT,0.0,1.0,ILOBOOL);
			y_jkt[j].push_back(var);
		}	
	}

	//ajout variables z_kt binaires
	vector<IloNumVarArray> z_kt;
	for(int k=0; k<d.cardM; ++k){
		IloNumVarArray var (env,d.cardT,0.0,1.0,ILOBOOL);
		z_kt.push_back(var);
	}

	// CONTRAINTES

	//contraintes sur ct
	/*IloNumArray bpt(env, d.nb_bp[0]-1);
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
	}*/
	for(int i=0; i<d.cardT; ++i){
		IloExpr pwct(env);
		for(int j=0; j<d.nb_bp[i];++j){
			pwct += d.pente[i][j]*xij[i][j]+(d.valbpt[i][j]-d.pente[i][j]*d.bpt[i][j])*pwd[i][j];
		}
		model.add(ct[i] == pwct);
	}

	// Pminj*pwdij <= xij
	for(int i=0; i<d.cardT; ++i){
		for(int j=0; j<d.nb_bp[i]; ++j){
			model.add(pwd[i][j]*d.bpt[i][j] <= xij[i][j] <= pwd[i][j]*d.bpt[i][j+1]-1);
		}
	}

	for(int i=0; i<d.cardT; ++i){
		model.add(IloSum(pwd[i]) == 1);
		model.add(IloSum(xij[i]) == xt[i]);
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
if(p.taches_avec_fenetre_temps == 1){
	// contrainte respect fenetre temps
	for(int j=0; j<d.cardJ; ++j){
        IloExpr s1(env);
		IloExpr s2(env);
        for(int k=0; k<d.cardM; ++k){
			for(int t=0; t<d.rj[j]; ++t){
				s1 += y_jkt[j][k][t];
			}
			for(int t=d.dj[j]; t<d.cardT; ++t){
				s2 += y_jkt[j][k][t];
			}
		}			
		model.add(s1 + s2 == 0);
	}
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

	// xt -st +st-1 - sum_j sum_k Djk*yjkt - sum_k Dk*zkt = 0
	
    
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
        if(t!=0) model.add(xt[t]-st[t]+st[t-1]-sum1-sum2 == 0);
        else model.add(xt[t]-st[t]+d.s0-sum1-sum2 == 0);
	}

    model.add(st[d.cardT-1]-d.s0 >= 0);



    auto start_time = chrono::steady_clock::now();
    cplex.solve();
	auto end_time = chrono::steady_clock::now();
	tps = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;



	/*for(int t=0; t<d.cardT; ++t){
		cout << "xt" << t << ": " <<cplex.getValue(xt[t]) << endl;
		cout << "ct" << t << ": " <<cplex.getValue(ct[t]) << endl;
		cout << "st" << t << ": " <<cplex.getValue(st[t]) << endl;
		for(int k=0; k<d.cardM; ++k){
			cout << "zkt" << k << t << ": " <<cplex.getValue(z_kt[k][t]) << endl;
			for(int j=0; j<d.cardJ; ++j){
			
			
				cout << "yjkt" << j << k<<t << ": " <<cplex.getValue(y_jkt[j][k][t]) << endl;
			}
		}
	
	}*/

	
    float sol;
	if(cplex.getStatus() == IloAlgorithm::Feasible){
		status = "Feasible";
		borneinf = cplex.getBestObjValue();
		sol = cplex.getObjValue();
		vector<float> variat;
		variat.push_back(roundd(cplex.getValue(st[0]),5) - d.s0);
		for(int t=1; t<d.cardT; ++t){
			variat.push_back(roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5));
		}
		relax = variat;
	}
    else if(cplex.getStatus() == IloAlgorithm::Optimal){
		status = "Optimal";
		borneinf = cplex.getBestObjValue();
		sol = cplex.getObjValue();
		vector<float> variat;
		variat.push_back(roundd(cplex.getValue(st[0]),5) - d.s0);
		for(int t=1; t<d.cardT; ++t){
			variat.push_back(roundd(cplex.getValue(st[t]) - cplex.getValue(st[t-1]),5));
			//cout << "var : "<<variat[t]<<endl;
		}
		relax = variat;
	}
	else{
		status = "Unknown";
		borneinf = cplex.getBestObjValue();
		sol = -1.0;
	}

    env.end();
    return sol;
}