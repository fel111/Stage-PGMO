#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "data_struct.h"
//#include "lotsizcontcom.h"
#include <ilcplex/ilocplex.h>
#include "lotsizcontCPX.h"

//#define SCIP_DEBUG
ILOSTLBEGIN

float lotsizcontcplex(data d){
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	cplex.setParam(IloCplex::Threads, 1);
	cplex.setParam(IloCplex::TiLim,300);

	//VARIABLES
	//ajout variables ct
	IloNumVarArray ct (env,d.cardT,0.0,IloInfinity,ILOFLOAT);
	
		
	//ajout variables st de stock
    IloNumVarArray st (env,d.cardT,0.0,d.Q,ILOFLOAT);


	//ajout variables zt binaires
	/*vector<vector<SCIP_VAR *> > ztj;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_VAR *> v;
		ztj.push_back(v);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_VAR * var;
			ztj[i].push_back(var);
			SCIPcreateVarBasic(scip, &ztj[i][j], ("ztj"+to_string(i)+to_string(j)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(scip,ztj[i][j]);		
		}	
	}*/

	//ajout variables rt de production, rt = sum_j xij
	IloNumVarArray xt (env,d.cardT,0.0,IloInfinity,ILOFLOAT);

	
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
	model.add(xt[0] - st[0] == d.dt[0]);
	for(int i=1; i<d.cardT; ++i){
         model.add(xt[i] + st[i-1] - st[i] == d.dt[i]);
	}


	//contrainte stmax - st0 >= 0
	model.add(st[d.cardT-1] >= d.s0);

		
    //IloCplex cplex(env);
    //cplex.extract(model);
    //cplex.exportModel("cplex_lotsizcont.lp");
    cplex.solve();

    float sol = cplex.getObjValue();

    env.end();
	return sol;

}