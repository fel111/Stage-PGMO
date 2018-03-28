#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "data_struct.h"
#include "compact.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
//int infini_int = numeric_limits<int>::max();


void modele_entier_compact(data d){

	// DONNEES

	/*//int d.cardT;
	int d.cardJ = 4;
	int d.cardM = 2;
	int d.cardR = 2;
	int d.s0 = 10;
	int d.cardT = 10; 
	int d.Q = 20;
	//vector<int> dt = {5,10,15,20,25,12,45,24,12,65};
	vector<vector<int> > d.bpt;
	vector<vector<float> > d.valbpt;
	vector<vector<float> > d.pente;
	vector<int> d.nb_bp; 
	for(int i=0;i<d.cardT;++i){
		vector<float> temp_d.pente = {3.0,2.0,1.0};
		vector<int> temp_d.bpt = {0,10,20,infini_int};
		vector<float> temp_d.valbpt = {0.0,30.0,50.0};
		d.pente.push_back(temp_d.pente);
		d.bpt.push_back(temp_d.bpt);
		d.valbpt.push_back(temp_d.valbpt);
		d.nb_bp.push_back(3);
	}

	vector<int> d.pj = {4,2,8,5};
	vector<vector<float> > d.cjr = {{20.0,0.4},{15.0,0.2},{35.0,1.5},{10.0,1.7}};
	vector<vector<float> > d.Ckr = {{300.0,10.5},{500.0,10.4}};
	vector<float> d.Dk = {10.0,15.0};
	vector<vector<float> > D = {{10.0,20.3},{12.1,30.2},{12.4,5.2},{17.1,54.6}};
*/


	// MODELE



	SCIP * scip;
	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPsetMessagehdlr(scip,NULL);
	SCIPcreateProb(scip, "Compact", NULL, NULL, NULL, NULL, NULL, NULL, NULL);


	//VARIABLES
	//ajout variables ct
	vector<SCIP_VAR *> ct;
	for(int i=0; i<d.cardT; ++i){
		SCIP_VAR * var;
		ct.push_back(var);
		SCIPcreateVarBasic(scip, &ct[i], ("ct"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 1.0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip,ct[i]);
	}

	//ajout variables xtij
	vector<vector<SCIP_VAR *> > xtij;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_VAR *> v;
		xtij.push_back(v);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_VAR * var;
			xtij[i].push_back(var);
			SCIPcreateVarBasic(scip, &(xtij[i][j]), ("xtij"+to_string(i)+to_string(j)).c_str(), 0, SCIP_REAL_MAX, 0.0, SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(scip,xtij[i][j]);		
		}	
	}

	/*//ajout variables xt temp
	vector<vector<SCIP_VAR *> > xtemp;
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
	for(int i=0; i<=d.cardT; ++i){
		SCIP_VAR * var;
		st.push_back(var);
		SCIPcreateVarBasic(scip, &st[i], ("st"+to_string(i)).c_str(), 0, d.Q, 0.0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip,st[i]);
	}



	//ajout variables pwd binaires
	vector<vector<SCIP_VAR *> > pwd;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_VAR *> v;
		pwd.push_back(v);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_VAR * var;
			pwd[i].push_back(var);
			SCIPcreateVarBasic(scip, &pwd[i][j], ("pwd"+to_string(i)+to_string(j)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(scip,pwd[i][j]);		
		}	
	}

	//ajout variables rt de production, rt = sum_j xij
	vector<SCIP_VAR *> xt;
	for(int i=0; i<d.cardT; ++i){
		SCIP_VAR * var;
		xt.push_back(var);
		SCIPcreateVarBasic(scip, &xt[i], ("xt"+to_string(i)).c_str(), 0, SCIP_REAL_MAX, 0, SCIP_VARTYPE_CONTINUOUS);
		SCIPaddVar(scip,xt[i]);
	}

	//ajout variables y_jkt binaires
	vector<vector<vector<SCIP_VAR *> > > y_jkt;
	for(int j=0; j<d.cardJ; ++j){
		vector<vector<SCIP_VAR *> > v;
		y_jkt.push_back(v);
		for(int k=0; k<d.cardM; ++k){
			vector<SCIP_VAR *> var;
			y_jkt[j].push_back(var);
			for(int t=0; t<d.cardT; ++t){
				SCIP_VAR * varr;
				y_jkt[j][k].push_back(varr);
				SCIPcreateVarBasic(scip, &y_jkt[j][k][t], ("y_jkt"+to_string(j)+to_string(k)+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
				SCIPaddVar(scip,y_jkt[j][k][t]);		
		
			}
		}	
	}

	//ajout variables z_kt binaires
	vector<vector<SCIP_VAR *> > z_kt;
	for(int k=0; k<d.cardM; ++k){
		vector<SCIP_VAR *> v;
		z_kt.push_back(v);
		for(int t=0; t<d.cardT; ++t){
			SCIP_VAR * var;
			z_kt[k].push_back(var);
			SCIPcreateVarBasic(scip, &z_kt[k][t], ("z_kt"+to_string(k)+to_string(t)).c_str(), 0, 1, 0.0, SCIP_VARTYPE_BINARY);
			SCIPaddVar(scip,z_kt[k][t]);		
		}	
	}

	// CONTRAINTES

	//contraintes sur ct
	vector<SCIP_CONS *> cons_ct;
	for(int i=0; i<d.cardT; ++i){
		SCIP_CONS * c;
		cons_ct.push_back(c);
		SCIPcreateConsLinear(scip, &cons_ct[i], ("cons_ct"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_ct[i], ct[i],1);		
		for(int j=0; j<d.nb_bp[i]; ++j){
			//SCIPaddCoefLinear(scip, cons_ct[i], xt[i][j],-d.pente[i][j]);
			SCIPaddCoefLinear(scip, cons_ct[i], xtij[i][j],-d.pente[i][j]);
			SCIPaddCoefLinear(scip, cons_ct[i], pwd[i][j],-d.valbpt[i][j]+d.pente[i][j]*d.bpt[i][j]);
		}
		SCIPaddCons(scip, cons_ct[i]);
	}

	// Pminj*pwdij <= xij
	vector<vector<SCIP_CONS *> > cons_pwd;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_pwd.push_back(c);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_pwd[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_pwd[i][j], ("cons_pwd_binf"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, 0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_pwd[i][j], xtij[i][j],1);
			SCIPaddCoefLinear(scip, cons_pwd[i][j], pwd[i][j],-d.bpt[i][j]);
			SCIPaddCons(scip, cons_pwd[i][j]);
		}
	}

	// xij < Pmaxj*Zij
	vector<vector<SCIP_CONS *> > cons_pwd2;
	for(int i=0; i<d.cardT; ++i){
		vector<SCIP_CONS *> c;
		cons_pwd2.push_back(c);
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIP_CONS * cons;
			cons_pwd2[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_pwd2[i][j], ("cons_pwd_bsup"+to_string(i)+to_string(j)).c_str(), 0, 0, 0, -SCIPinfinity(scip), 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_pwd2[i][j], xtij[i][j],1);
			SCIPaddCoefLinear(scip, cons_pwd2[i][j], pwd[i][j],-d.bpt[i][j+1]+1);
			SCIPaddCons(scip, cons_pwd2[i][j]);
		}
	}	

	// sum pwdij = 1
	vector<SCIP_CONS *> cons_pwdsum;
	for(int i=0; i<d.cardT; ++i){
		SCIP_CONS * c;
		cons_pwdsum.push_back(c);
		SCIPcreateConsLinear(scip, &cons_pwdsum[i], ("cons_pwdsum"+to_string(i)).c_str(), 0, 0, 0, 1, 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_pwdsum[i], pwd[i][j],1);
		}
		SCIPaddCons(scip, cons_pwdsum[i]);
	}
	
	// sum xtij = xt
	vector<SCIP_CONS *> cons_xt;
	for(int i=0; i<d.cardT; ++i){
		SCIP_CONS * c;
		cons_xt.push_back(c);
		SCIPcreateConsLinear(scip, &cons_xt[i], ("cons_xt"+to_string(i)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_xt[i], xt[i],1);				
		for(int j=0; j<d.nb_bp[i]; ++j){
			SCIPaddCoefLinear(scip, cons_xt[i], xtij[i][j],-1);
		}
		SCIPaddCons(scip, cons_xt[i]);
	}
	
	// sum y_jkt <= 1
	vector<vector<SCIP_CONS *> > cons_y_jkt;
	for(int j=0; j<d.cardJ; ++j){
		vector<SCIP_CONS *> c;
		cons_y_jkt.push_back(c);
		for(int t=0; t<d.cardT; ++t){
			SCIP_CONS * cons;
			cons_y_jkt[j].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_y_jkt[j][t], ("cons_y_jkt"+to_string(j)+to_string(t)).c_str(), 0, 0, 0, -SCIPinfinity(scip), 1, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
			for(int k=0; k<d.cardM; ++k){
				SCIPaddCoefLinear(scip, cons_y_jkt[j][t], y_jkt[j][k][t],1);
			}			
			SCIPaddCons(scip, cons_y_jkt[j][t]);
		}
	}

	// sum_k sum_t y_jkt = d.pj
	vector<SCIP_CONS *> cons_y_jkt2;
	for(int j=0; j<d.cardJ; ++j){
		SCIP_CONS * c;
		cons_y_jkt2.push_back(c);
		SCIPcreateConsLinear(scip, &cons_y_jkt2[j], ("cons_y_jkt2"+to_string(j)).c_str(), 0, 0, 0, d.pj[j], d.pj[j], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		for(int k=0; k<d.cardM; ++k){
			for(int t=0; t<d.cardT; ++t){
				SCIPaddCoefLinear(scip, cons_y_jkt2[j], y_jkt[j][k][t],1);
			}
		}			
		SCIPaddCons(scip, cons_y_jkt2[j]);
	}

	// sum_j d.cjr*yjkt <= d.Ckr
	vector<vector<vector<SCIP_CONS *> > > cons_ckr;
	for(int k=0; k<d.cardM; ++k){
		vector<vector<SCIP_CONS *> > c;
		cons_ckr.push_back(c);
		for(int r=0; r<d.cardR; ++r){
			vector<SCIP_CONS *> cons;
			cons_ckr[k].push_back(cons);
			for(int t=0; t<d.cardT; ++t){
				SCIP_CONS * co;
				cons_ckr[k][r].push_back(co);
				SCIPcreateConsLinear(scip, &cons_ckr[k][r][t], ("cons_ckr"+to_string(k)+to_string(r)+to_string(t)).c_str(), 0, 0, 0, -SCIPinfinity(scip), d.Ckr[k][r], TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
				for(int j=0; j<d.cardJ; ++j){
					SCIPaddCoefLinear(scip, cons_ckr[k][r][t], y_jkt[j][k][t],d.cjr[j][r]);
				}
				SCIPaddCons(scip, cons_ckr[k][r][t]);			
			}
		}
	}

	// zkt - yjkt >= 0
	vector<vector<vector<SCIP_CONS *> > > cons_zkt;
	for(int k=0; k<d.cardM; ++k){
		vector<vector<SCIP_CONS *> > c;
		cons_zkt.push_back(c);
		for(int j=0; j<d.cardJ; ++j){
			vector<SCIP_CONS *> cons;
			cons_zkt[k].push_back(cons);
			for(int t=0; t<d.cardT; ++t){
				SCIP_CONS * co;
				cons_zkt[k][j].push_back(co);
				SCIPcreateConsLinear(scip, &cons_zkt[k][j][t], ("cons_zkt"+to_string(k)+to_string(j)+to_string(t)).c_str(), 0, 0, 0, 0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
				SCIPaddCoefLinear(scip, cons_zkt[k][j][t], y_jkt[j][k][t],-1);
				SCIPaddCoefLinear(scip, cons_zkt[k][j][t], z_kt[k][t],1);
				SCIPaddCons(scip, cons_zkt[k][j][t]);			
			}
		}
	}

	// xt + st - st-1 - sum_j sum_k Djk*yjkt - sum_k d.Dk*zkt = 0
	vector<SCIP_CONS *> cons_equil;
	SCIP_CONS * c;
	cons_equil.push_back(c);
	SCIPcreateConsLinear(scip, &cons_equil[0], "cons_equil0", 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(scip, cons_equil[0], xt[0],1);
	SCIPaddCoefLinear(scip, cons_equil[0], st[0],1);		
	for(int j=0; j<d.cardJ; ++j){
		for(int k=0; k<d.cardM; ++k){
			SCIPaddCoefLinear(scip, cons_equil[0], y_jkt[j][k][0],-d.Djk[j][k]);
		}
	}
	for(int k=0; k<d.cardM; ++k){
		SCIPaddCoefLinear(scip, cons_equil[0], z_kt[k][0],-d.Dk[k]);
	}			
	SCIPaddCons(scip, cons_equil[0]);
	for(int t=1; t<d.cardT; ++t){
		SCIP_CONS * c;
		cons_equil.push_back(c);
		SCIPcreateConsLinear(scip, &cons_equil[t], ("cons_equil"+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
		SCIPaddCoefLinear(scip, cons_equil[t], xt[t],1);
		SCIPaddCoefLinear(scip, cons_equil[t], st[t],1);
		SCIPaddCoefLinear(scip, cons_equil[t], st[t-1],-1);		
		for(int j=0; j<d.cardJ; ++j){
			for(int k=0; k<d.cardM; ++k){
				SCIPaddCoefLinear(scip, cons_equil[t], y_jkt[j][k][t],-d.Djk[j][k]);
			}
		}
		for(int k=0; k<d.cardM; ++k){
			SCIPaddCoefLinear(scip, cons_equil[t], z_kt[k][t],-d.Dk[k]);
		}			
		SCIPaddCons(scip, cons_equil[t]);
	}

	//contrainte stmax - st0 >= 0
	SCIP_CONS * cons_st;
	SCIPcreateConsLinear(scip, &cons_st, "cons_st", 0, 0, 0, 0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(scip, cons_st, st[d.cardT-1], 1);
	SCIPaddCoefLinear(scip, cons_st, st[0], -1);
	SCIPaddCons(scip, cons_st);

	//contrainte st0 = d.s0
	SCIP_CONS * cons_st0;
	SCIPcreateConsLinear(scip, &cons_st0, "cons_st0", 0, 0, 0, d.s0, d.s0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE);  
	SCIPaddCoefLinear(scip, cons_st0, st[0], 1);
	SCIPaddCons(scip, cons_st0);

	//FILE * filed;
	//filed = fopen("sol_comp", "w");
	//SCIPprintOrigProblem(scip, filed, "lp", false);
	

	SCIPsolve(scip);
	//RECUPERATION
	SCIP_SOL * sol = SCIPgetBestSol(scip);
	cout << "obj = " << SCIPgetSolOrigObj(scip,sol) << endl << endl;
	//SCIPprintBestSol(scip,filed,true); 	
	//fclose(filed);
	
}
