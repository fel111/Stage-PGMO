#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <ilcplex/cplex.h>
using namespace std;

float infini = numeric_limits<float>::infinity();


float p_t(int t,vector<int> nb_bp,vector<vector <int> > bpt,vector<vector <int> > valbpt,vector<vector <float> > pente, float x){
	if (x==0.0) return 0.0;	
	for (int i=0; i<nb_bp[t]; ++i){
		if((x>bpt[t][i])&&(x<=bpt[t][i+1])) 
			return pente[t][i]*x+valbpt[t][i];
	}
}


void lotSizComp(int cardT, int qmax, vector<int> dt, vector<vector<int> > bpt, vector<vector<int> > valbpt, vector<vector <float> > pente,vector<int> nb_bp){
	SCIP * scip;
	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPcreateProb(scip, "LotSizComp", 0, 0, 0, 0, 0, 0, 0);


	//VARIABLES
	//ajout variables ct
	vector<SCIP_VAR *> ct;
	for(int i=0; i<cardT; ++i){
		SCIP_VAR * var;
		ct.push_back(var);
		SCIPcreateVarBasic(scip, &ct[i], ("ct"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 1.0, SCIP_VARTYPE_CONTINUOUS);
	}

	//ajout variables rt
	vector<SCIP_VAR *> xt;
	for(int i=0; i<cardT; ++i){
		SCIP_VAR * var;
		xt.push_back(var);
		SCIPcreateVarBasic(scip, &xt[i], ("xt"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_CONTINUOUS);
	}
	
	//ajout variables st de stock
	vector<SCIP_VAR *> st;
	for(int i=0; i<cardT; ++i){
		SCIP_VAR * var;
		st.push_back(var);
		SCIPcreateVarBasic(scip, &st[i], ("st"+to_string(i)).c_str(), 0, qmax, 1.0, SCIP_VARTYPE_CONTINUOUS);
	}


	//CONTRAINTES
	//contraintes sur ct
	vector<SCIP_CONS *> cons_ct;
	for(int i=0; i<cardT; ++i){
		SCIP_CONS * cons;
		cons_ct.push_back(cons);
		SCIPcreateConsLinear(scip, &cons_ct[i], ("cons_ct"+to_string(i)).c_str(), 0, 0, 0, p_t(i,nb_bp,bpt,valbpt,pente,xt[i]), p_t(i,nb_bp,bpt,valbpt,pente,xt[i]), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_ct[i], ct[i],1);
		SCIPaddCons(scip, cons_ct[i]);
	}
	
	//contraintes bilan energetique
	/*vector<SCIP_CONS *> cons_dt;
	for(int i=1; i<cardT; ++i){
		SCIP_CONS * cons;
		cons_dt.push_back(cons);
		SCIPcreateConsLinear(scip, &cons_dt[i], ("cons_dt"+to_string(i)).c_str(), 0, 0, 0, dt[i], dt[i], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_dt[i], st[i],1);
		SCIPaddCoefLinear(scip, cons_dt[i], st[i-1],-1);
		SCIPaddCoefLinear(scip, cons_dt[i], rt[i],1);		
		SCIPaddCons(scip, cons_dt[i]);
	}*/

	//contrainte stmax = 0
	/*SCIP_CONS * cons_stmax;
	SCIPcreateConsLinear(scip,&cons_stmax,"cons_qtmax",0,0,0,0,0,TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  	
	SCIPaddCoefLinear(scip,cons_stmax,st[cardT-1],1);
	SCIPaddCons(scip,cons_stmax);*/

	//contrainte evolution stock
	/*vector<SCIP_CONS *> cons_st;
	for(int i=0; i<cardT-1; ++i){
		SCIP_CONS * cons;
		cons_st.push_back(cons);
		SCIPcreateConsLinear(scip, &cons_st[i], ("cons_st"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_st[i], st[i+1],1);
		SCIPaddCoefLinear(scip, cons_st[i], st[i],-1);
		SCIPaddCoefLinear(scip, cons_st[i], rt[i],1);		
		SCIPaddCons(scip, cons_st[i]);
	}*/


	//RESOLUTION
	SCIPsolve(scip);
	
	//RECUPERATION
	SCIP_SOL * sol = SCIPgetBestSol(scip);
	cout << "obj = " << SCIPgetSolOrigObj(scip,sol) << endl << endl;
	/*cout << "etat de la charge :" << endl;	
	for(int i=0; i<cardT; ++i){
		cout << "s"<<i<<" = "<<SCIPgetSolVal(scip,sol,st[i]) << endl;
	}
	cout << endl;
	cout << "Cout energie : "<<endl;
	for(int i=0; i<cardT; ++i){
		cout << "ct"<<i<<" = "<<SCIPgetSolVal(scip,sol,ct[i]) << endl;
	}*/
	

}

int main(){
	int cardT = 5; 
	int qmax = 100;
	vector<int> dt = {5,10,15,20,25};
	vector<vector<int> > bpt;
	vector<vector<int> > valbpt;
	vector<vector<float> > pente;
	vector<int> nb_bp = {3,3,3,3,3};
	for(int i=0;i<cardT;++i){
		vector<float> temp_pente = {1.0,2.0,infini};
		vector<int> temp_bpt = {0,100,200};
		vector<int> temp_valbpt = {0,100,300};
		pente.push_back(temp_pente);
		bpt.push_back(temp_bpt);
		valbpt.push_back(temp_valbpt);
	}

	lotSizComp(cardT,qmax,dt,bpt,valbpt,pente,nb_bp);

	return 0;
}

