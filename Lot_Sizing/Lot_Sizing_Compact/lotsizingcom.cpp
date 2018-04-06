#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <chrono>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "data_struct.h"
#include "lotsizingcom.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;

//int infini_int = numeric_limits<int>::max();


//a inclure dans un header ensuite
/*void lecture_pwd(string file, vector<vector<int> >& d.bpt, vector<vector<float> >& d.valbpt, vector<vector<float> >& d.pente){

	//char temp;
	//vector< vector<int> > conflitstemp;
	
	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		int bp;
		float valbp;
		float p;
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> bp >> valbp >> p;
			d.bpt[t].push_back(bp);
			d.valbpt[t].push_back(valbp);
			d.pente[t].push_back(p);
			//cout << "OK" << endl;
		}
		fichier.close();
		//cout << " done "<<endl;
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}

}

//a inclure dans un header ensuite
void lecture_demande(string file, vector<int> & dt){

	//char temp;
	//vector< vector<int> > conflitstemp;
	
	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		int d;
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> d;
			dt.push_back(d);
			//cout << "dt["<<t<<"] : "<<dt[t]<<endl;
		}
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}*/


void lotsizcomp(data d, int choix){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();

	SCIP * scip;
	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPsetMessagehdlr(scip,NULL);
	SCIPcreateProb(scip, "LotSizingComp", NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	vector<int> dt = dtToInt(d,choix);

	//VARIABLES
	//ajout variables ct
	vector<SCIP_VAR *> ct;
	for(int i=0; i<d.cardT; ++i){
		SCIP_VAR * var;
		ct.push_back(var);
		SCIPcreateVarBasic(scip, &ct[i], ("ct"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 1.0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip,ct[i]);
	}

	//ajout variables xt
	vector<vector<SCIP_VAR *> > xt;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_VAR *> v;
		xt.push_back(v);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_VAR * var;
			xt[i].push_back(var);
			SCIPcreateVarBasic(scip, &(xt[i][j]), ("xt"+to_string(i)+to_string(j)).c_str(), 0, SCIP_REAL_MAX, 0.0, SCIP_VARTYPE_INTEGER);
			SCIPaddVar(scip,xt[i][j]);		
		}	
	}

	//ajout variables xt temp
	/*vector<vector<SCIP_VAR *> > xtemp;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_VAR *> v;
		xtemp.push_back(v);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_VAR * var;
			xtemp[i].push_back(var);
			SCIPcreateVarBasic(scip, &(xtemp[i][j]), ("xtemp"+to_string(i)+to_string(j)).c_str(), 0, SCIP_REAL_MAX, 0.0, SCIP_VARTYPE_INTEGER);
			SCIPaddVar(scip,xtemp[i][j]);		
		}	
	}*/

	
	//ajout variables st de stock
	vector<SCIP_VAR *> st;
	for(int i=0; i<d.cardT; ++i){
		SCIP_VAR * var;
		st.push_back(var);
		SCIPcreateVarBasic(scip, &st[i], ("st"+to_string(i)).c_str(), 0, d.Q, 0.0, SCIP_VARTYPE_INTEGER);
		SCIPaddVar(scip,st[i]);
	}



	//ajout variables zt binaires
	vector<vector<SCIP_VAR *> > ztj;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_VAR *> v;
		ztj.push_back(v);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_VAR * var;
			ztj[i].push_back(var);
			SCIPcreateVarBasic(scip, &ztj[i][j], ("ztj"+to_string(i)+to_string(j)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(scip,ztj[i][j]);		
		}	
	}

	//ajout variables rt de production, rt = sum_j xij
	vector<SCIP_VAR *> rt;
	for(int i=0; i<d.cardT; ++i){
		SCIP_VAR * var;
		rt.push_back(var);
		SCIPcreateVarBasic(scip, &rt[i], ("rt"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_INTEGER);
		SCIPaddVar(scip,rt[i]);
	}

	
	//CONTRAINTES
	//contraintes sur ct
	vector<SCIP_CONS *> cons_ct;
	for(int i=0; i<d.cardT; ++i){
		SCIP_CONS * c;
		cons_ct.push_back(c);
		SCIPcreateConsLinear(scip, &cons_ct[i], ("cons_ct"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_ct[i], ct[i],1);		
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_ct[i], xt[i][j],-d.pente[i][j]);
			//SCIPaddCoefLinear(scip, cons_ct[i], xtemp[i][j],-d.pente[i][j]);
			SCIPaddCoefLinear(scip, cons_ct[i], ztj[i][j],-d.valbpt[i][j]+d.pente[i][j]*d.bpt[i][j]);
		}
		SCIPaddCons(scip, cons_ct[i]);
	}

	// Pminj*Zij <= xij
	vector<vector<SCIP_CONS *> > cons_z;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_z.push_back(c);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_z[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_z[i][j], ("cons_z_binf"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, 0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_z[i][j], xt[i][j],1);
			SCIPaddCoefLinear(scip, cons_z[i][j], ztj[i][j],-d.bpt[i][j]);
			SCIPaddCons(scip, cons_z[i][j]);
		}
	}

	// xij < Pmaxj*Zij
	vector<vector<SCIP_CONS *> > cons_z2;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_z2.push_back(c);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_z2[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_z2[i][j], ("cons_z_bsup"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, -SCIPinfinity(scip), 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_z2[i][j], xt[i][j],1);
			SCIPaddCoefLinear(scip, cons_z2[i][j], ztj[i][j],-d.bpt[i][j+1]+1);
			SCIPaddCons(scip, cons_z2[i][j]);
		}
	}

	// sum zij = 1
	vector<SCIP_CONS *> cons_zsum;
	for(int i=0; i<d.cardT; ++i){
		SCIP_CONS * c;
		cons_zsum.push_back(c);
		SCIPcreateConsLinear(scip, &cons_zsum[i], ("cons_zsum"+to_string(i)).c_str(), 0, 0, 0, 1, 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_zsum[i], ztj[i][j],1);
		}
		SCIPaddCons(scip, cons_zsum[i]);
	}

	// sum xt_ij = rt
	vector<SCIP_CONS *> cons_rt;
	for(int i=0; i<d.cardT; ++i){
		SCIP_CONS * c;
		cons_rt.push_back(c);
		SCIPcreateConsLinear(scip, &cons_rt[i], ("cons_rt"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_rt[i], rt[i],1);				
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_rt[i], xt[i][j],-1);
		}
		SCIPaddCons(scip, cons_rt[i]);
	}

	// contrainte xtemp_ij = x_ij - d.bpt_ij
	/*vector<vector<SCIP_CONS *> > cons_xtemp;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_xtemp.push_back(c);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_xtemp[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_xtemp[i][j], ("cons_xtemp"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_xtemp[i][j], xt[i][j],1);
			SCIPaddCoefLinear(scip, cons_xtemp[i][j], xt[i][j],-1);
			SCIPaddCoefLinear(scip, cons_xtemp[i][j], ztj[i][j],d.bpt[i][j]);
			SCIPaddCons(scip, cons_xtemp[i][j]);
		}
	}*/

	
	//contraintes bilan energetique
	vector<SCIP_CONS *> cons_dt;
	SCIP_CONS * cons;
	cons_dt.push_back(cons);
	SCIPcreateConsLinear(scip, &cons_dt[0], "cons_dt0", 0, 0, 0, dt[0], dt[0], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(scip, cons_dt[0], st[0],-1);
	SCIPaddCoefLinear(scip, cons_dt[0], rt[0],1);		
	SCIPaddCons(scip, cons_dt[0]);
	for(int i=1; i<d.cardT; ++i){
		SCIP_CONS * cons;
		cons_dt.push_back(cons);
		SCIPcreateConsLinear(scip, &cons_dt[i], ("cons_dt"+to_string(i)).c_str(), 0, 0, 0, dt[i], dt[i], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_dt[i], st[i],-1);
		SCIPaddCoefLinear(scip, cons_dt[i], st[i-1],1);
		SCIPaddCoefLinear(scip, cons_dt[i], rt[i],1);		
		SCIPaddCons(scip, cons_dt[i]);
	}


	//contrainte stmax = 0
	SCIP_CONS * cons_stmax;
	SCIPcreateConsLinear(scip, &cons_stmax, "cons_stmax", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(scip, cons_stmax, st[d.cardT-1], 1);
	SCIPaddCons(scip, cons_stmax);
	
	//contrainte st0 = d.s0
	SCIP_CONS * cons_st0;
	SCIPcreateConsLinear(scip, &cons_st0, "cons_st0", 0, 0, 0, d.s0, d.s0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(scip, cons_st0, st[0], 1);
	SCIPaddCons(scip, cons_st0);

	//RESOLUTION
	//SCIP_CALL(status = SCIPsolve(scip));
	SCIPsolve(scip);
	//RECUPERATION
	SCIP_SOL * sol = SCIPgetBestSol(scip);
	cout << "obj comp = " << SCIPgetSolOrigObj(scip,sol) << endl;

	/*FILE * filed;
	filed = fopen("SolLotSizComp", "a+");
	SCIPprintOrigProblem(scip, filed, "lp", false);
	fclose(filed);*/



	/*cout << "production :" << endl;	
	for(int i=0; i<d.cardT; ++i){
		cout << "x"<<i<<" = "<<SCIPgetSolVal(scip,sol,rt[i]) << endl;
	}
	cout << endl;
	cout << "Cout energie : "<<endl;
	for(int i=0; i<d.cardT; ++i){
		cout << "ct"<<i<<" = "<<SCIPgetSolVal(scip,sol,ct[i]) << endl;
	}
	cout << "stock :" << endl;	
	for(int i=0; i<d.cardT; ++i){
		cout << "s"<<i<<" = "<<SCIPgetSolVal(scip,sol,st[i]) << endl;
	}
	cout << endl;

	for(int i=0;i<d.cardT;++i){
		for(int j=0; j<d.nb_bp[i]; ++j){
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



/*
int main(){
	int cardT = 10; 
	int d.Q = 20;
	vector<int> dt = {5,10,15,20,25,12,45,24,12,65};
	vector<vector<int> > d.bpt;
	vector<vector<float> > d.valbpt;
	vector<vector<float> > d.pente;
	vector<int> d.nb_bp = {3,3,3,3,3,3,3,3,3,3};
	for(int i=0;i<cardT;++i){
		vector<float> temp_d.pente = {3.0,2.0,1.0};
		vector<int> temp_d.bpt = {0,10,20,infini_int};
		vector<float> temp_d.valbpt = {0.0,30.0,50.0};
		d.pente.push_back(temp_d.pente);
		d.bpt.push_back(temp_d.bpt);
		d.valbpt.push_back(temp_d.valbpt);
	}
	int cardT = 100;
	int d.Q = 50;
	vector<int> d.nb_bp (cardT, 10);
	vector<vector<int> > d.bpt;
	vector<vector<float> > d.valbpt;
	vector<vector<float> > d.pente;
	for(int i=0; i<cardT; ++i){
		vector<int> bp;
		vector<float> valbp;
		vector<float> p;
		d.bpt.push_back(bp);
		d.valbpt.push_back(valbp);
		d.pente.push_back(p);	
	}
	vector<int> dt;
	
	lecture_demande("demande_100_200", dt);
	lecture_pwd("pwd_100_10_3_100", d.bpt, d.valbpt, d.pente);

	for(int k=0; k<cardT;++k){
		d.bpt[k].push_back(infini_int);
	}

	lotSizComp(cardT,d.Q,dt,d.bpt,d.valbpt,d.pente,d.nb_bp);

	return 0;
}*/

