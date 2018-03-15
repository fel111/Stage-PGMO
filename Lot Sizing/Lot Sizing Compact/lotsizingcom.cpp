#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <ilcplex/cplex.h>
#define SCIP_DEBUG
using namespace std;

float infini = numeric_limits<float>::infinity();


/*float p_t(int t,vector<int> nb_bp,vector<vector <int> > bpt,vector<vector <int> > valbpt,vector<vector <float> > pente, float x){
	if (x==0.0) return 0.0;	
	for (int i=0; i<nb_bp[t]; ++i){
		if((x>bpt[t][i])&&(x<=bpt[t][i+1])) 
			return pente[t][i]*x+valbpt[t][i];
	}
}*/


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

	//ajout variables xt
	vector<vector<SCIP_VAR *> > xt;
	for(int i=0; i<cardT; ++i){
		vector<SCIP_VAR *> v;
		xt.push_back(v);
		for(int j=0; j<nb_bp[i]; ++j){
			SCIP_VAR * var;
			xt[i].push_back(var);
			SCIPcreateVarBasic(scip, &xt[i][j], ("xt"+to_string(i)+to_string(j)).c_str(), 0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_INTEGER);
		}	
	}
	
	//ajout variables st de stock
	vector<SCIP_VAR *> st;
	for(int i=0; i<cardT; ++i){
		SCIP_VAR * var;
		st.push_back(var);
		SCIPcreateVarBasic(scip, &st[i], ("st"+to_string(i)).c_str(), 0, qmax, 0, SCIP_VARTYPE_INTEGER);
	}

	/*SCIP_VAR ** st = (SCIP_VAR**) malloc(cardT*sizeof(SCIP_VAR*));
	for(int i=0; i<cardT; ++i){
		SCIP_VAR * var;
		st[i] = (var);
		SCIPcreateVarBasic(scip, &st[i], ("st"+to_string(i)).c_str(), 0, qmax, 0, SCIP_VARTYPE_INTEGER);
	}*/

	//ajout variables zt binaires
	vector<vector<SCIP_VAR *> > ztj;
	for(int i=0; i<cardT; ++i){
		vector<SCIP_VAR *> v;
		ztj.push_back(v);
		for(int j=0; j<nb_bp[i]; ++j){
			SCIP_VAR * var;
			ztj[i].push_back(var);
			SCIPcreateVarBasic(scip, &ztj[i][j], ("ztj"+to_string(i)+to_string(j)).c_str(), 0, 1, 0, SCIP_VARTYPE_BINARY);
		}	
	}

	//ajout variables rt de production, rt = sum_j xij
	vector<SCIP_VAR *> rt;
	for(int i=0; i<cardT; ++i){
		SCIP_VAR * var;
		rt.push_back(var);
		SCIPcreateVarBasic(scip, &rt[i], ("rt"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_INTEGER);
	}


	//CONTRAINTES
	//contraintes sur ct
	/*vector<SCIP_CONS *> cons_ct;
	for(int i=0; i<cardT; ++i){
		SCIP_CONS * c;
		cons_ct.push_back(c);
		SCIPcreateConsLinear(scip, &cons_ct[i], ("cons_ct"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_ct[i], ct[i],1);		
		for(int j=0; j<nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_ct[i], xt[i][j],-pente[i][j]);
			SCIPaddCoefLinear(scip, cons_ct[i], ztj[i][j],-valbpt[i][j]);
		}
		SCIPaddCons(scip, cons_ct[i]);
	}

	// Pminj*Zij <= xij
	vector<vector<SCIP_CONS *> > cons_z;
	for(int i=0; i<cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_z.push_back(c);
		for(int j=0; j<nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_z[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_z[i][j], ("cons_z_binf"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, 0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_z[i][j], xt[i][j],1);
			SCIPaddCoefLinear(scip, cons_z[i][j], ztj[i][j],-bpt[i][j]);
			SCIPaddCons(scip, cons_z[i][j]);
		}
	}

	// xij < Pmaxj*Zij
	vector<vector<SCIP_CONS *> > cons_z2;
	for(int i=0; i<cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_z2.push_back(c);
		for(int j=0; j<nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_z2[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_z2[i][j], ("cons_z_bsup"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, -SCIPinfinity(scip), 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_z2[i][j], xt[i][j],1);
			SCIPaddCoefLinear(scip, cons_z2[i][j], ztj[i][j],-bpt[i][j]+1);
			SCIPaddCons(scip, cons_z2[i][j]);
		}
	}

	// sum zij = 1
	vector<SCIP_CONS *> cons_zsum;
	for(int i=0; i<cardT; ++i){
		SCIP_CONS * c;
		cons_zsum.push_back(c);
		SCIPcreateConsLinear(scip, &cons_zsum[i], ("cons_zsum"+to_string(i)).c_str(), 0, 0, 0, 1, 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		for(int j=0; j<nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_zsum[i], ztj[i][j],1);
		}
		SCIPaddCons(scip, cons_zsum[i]);
	}

	// sum xt_ij = rt
	vector<SCIP_CONS *> cons_rt;
	for(int i=0; i<cardT; ++i){
		SCIP_CONS * c;
		cons_rt.push_back(c);
		SCIPcreateConsLinear(scip, &cons_rt[i], ("cons_rt"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_rt[i], rt[i],1);				
		for(int j=0; j<nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_rt[i], xt[i][j],-1);
		}
		SCIPaddCons(scip, cons_rt[i]);
	}

	
	//contraintes bilan energetique
	/*vector<SCIP_CONS *> cons_dt;
	for(int i=1; i<cardT; ++i){
		SCIP_CONS * cons;
		cons_dt.push_back(cons);
		SCIPcreateConsLinear(scip, &cons_dt[i-1], ("cons_dt"+to_string(i)).c_str(), 0, 0, 0, dt[i], dt[i], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_dt[i-1], st[i],-1);
		SCIPaddCoefLinear(scip, cons_dt[i-1], st[i-1],1);
		SCIPaddCoefLinear(scip, cons_dt[i-1], rt[i],1);		
		SCIPaddCons(scip, cons_dt[i-1]);
	}

	//contrainte stmax = 0
	SCIP_CONS * cons_stmax;
	SCIPcreateConsLinear(scip,&cons_stmax,"cons_qtmax",0,NULL,NULL,0,0,TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  	
	SCIPaddCoefLinear(scip,cons_stmax,st[cardT-1],1);
	SCIPaddCons(scip,cons_stmax);*/

	
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

	/*for(int i=0;i<cardT;++i){
		for(int j=0; j<nb_bp[i]; ++j){
			SCIPreleaseVar(&xt[i][j]);
			SCIPreleaseVar(&ztj[i][j]);
		}
		SCIPreleaseVar(&rt[i]);
		SCIPreleaseVar(&st[i]);
		SCIPreleaseVar(&ct[i]);
	}*/
	//SCIPreleaseVar(&xt);
	SCIPfree(&scip);
	

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

