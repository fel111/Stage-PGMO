#include <iostream>
#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
//int infini_int = numeric_limits<int>::max();


struct data{
	// PWD
	vector<int> nb_bp;
	vector<vector <int> > bpt;
	vector< vector<float> > valbpt;
	vector<vector <float> > pente;

	// Data Ordo
	int cardT, cardJ, cardM, cardR, s0, Q;
	vector<int> pj;
	vector<vector<float> > cjr;
	vector<vector<float> > Ckr;
	vector<float> Dk;
	vector<vector<float> > Djk;
	vector<vector<int> > ri;
	vector<vector<int> > di;
};

data init(){
	data d;
	//int cardT;
	d.cardJ = 4;
	d.cardM = 2;
	d.cardR = 2;
	d.s0 = 10;
	d.cardT = 10; 
	d.Q = 20;
	//vector<int> dt = {5,10,15,20,25,12,45,24,12,65};
	for(int i=0;i<d.cardT;++i){
		vector<float> temp_pente = {3.0,2.0,1.0};
		vector<int> temp_bpt = {0,10,20,infini_int};
		vector<float> temp_valbpt = {0.0,30.0,50.0};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(3);
	}

	d.pj = {4,2,8,5};
	d.cjr = {{20.0,0.4},{15.0,0.2},{35.0,1.5},{10.0,1.7}};
	d.Ckr = {{300.0,10.5},{500.0,10.4}};
	d.Dk = {10.0,15.0};
	d.Djk = {{10.0,20.3},{12.1,30.2},{12.4,5.2},{17.1,54.6}};

	return d;
}





float p_t(int t, float x, data d){
	//cout << "x: "<<x<<endl;	
	//if (x==0) return 0.0;	
	for (int i=0; i<d.nb_bp[t]; ++i){
		//cout<<bpt[t][i]<<"<"<<x<<"<="<<bpt[t][i+1]<<endl;
		if((x>=d.bpt[t][i])&&(x<d.bpt[t][i+1])){
			//cout << bpt[t][i] << endl;
			//cout << "valreel" <<  pente[t][i]*x+valbpt[t][i]<<endl;
			return d.pente[t][i]*(x-d.bpt[t][i])+d.valbpt[t][i];
			//return pente[t][i]*x+valbpt[t][i];
		}	
	}
	return 10000000.0;
}


/*vector<vector<int> > firstSol(data d){
	<vector<float> cl (d.cardT, 0.0);
	vector<vector<int> > res (d.cardT, ;
	//vector<vector<vector<int> > > yjkt (d.cardJ, vector<vector<int> > (d.cardM, vector<int> (d.cardT, 0)));
	
	int nbtaches = 0;	
	for(int t=0; t<d.cardT; ++t){
		for(int j=0; j<d.cardJ; ++j){
			for(int k=0; k<d.cardM; ++k){
				bool test = true;	
				for(int r=0; r<d.cardR; ++r){
					if((d.cjr[j][r]<=ress[k][r])&&(d.ri[j]>=t)&&(t+d.pj[j]<=d.dj[j])){
						
						cl[fixed] += p_t(t,d.Djk[j][k],d);
						yjkt[j][k][t] = 1;
					}else test = test && false;	
				}
				if(test){
					for(int r=0; r<d.cardR; ++r) ress[k][r] = ress[k][r]-d.cjr[j][r];
					
				}
			}
		}
		
	}
	
	
}*/

vector<int> firstSol(data d){
	vector<int> res;
	for(int j=0; j<d.cardJ; ++j){
		for(int t=0; t<d.cardT; ++t){
			for(int k=0; k<d.cardM; ++k){
				bool test = true;
				for(int r=0; r<d.cardR; ++r){
					if(!((d.cjr[j][r]<=d.Ckr[k][r])&&(d.ri[j]>=t)&&(t+d.pj[j]<d.dj[j]))){
						test = test && false;
					}	
				}
				if(test){
					res.push_back(j);
					return res;
					
				}
			}
		}
		
	}
	return res.push_back(-1);

}


int main(){
	
	
	data d = init();
	// MODELE
	vector<int> firstJ = firstSol(d);
	vector<vector<int> > Al = {firstJ};
	<vector<float> cl;
	(for int i=0; i<Al.size();++i){
		float cout = 0.0;
		for(int j=0; j<Al[i].size(); ++j){
			cout += p_t(0,Al[i][j],d);
		}
		cl.push_back(cout);
	}
	int L = 1;
	


	SCIP * scip;
	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPcreateProb(scip, "PMR", NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	
	// VARIABLES

	//ajout variables y_lt
	vector<vector<SCIP_VAR *> > y_lt;
	for(int l=0; l<L; ++l){
		vector<SCIP_VAR *> v;
		y_lt.push_back(v);
		for(int t=0; t<d.cardT; ++t){
			SCIP_VAR * var;
			y_lt[l].push_back(var);
			SCIPcreateVarBasic(scip, &(y_lt[l][t]), ("y"+to_string(l)+to_string(t)).c_str(), 0, SCIP_REAL_MAX, cl[l], SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(scip,y_lt[l][t]);		
		}	
	}

	//ajout variables x_it
	vector<vector<SCIP_VAR *> > x_it;
	for(int i=0; i<d.cardJ; ++i){
		vector<SCIP_VAR *> v;
		x_it.push_back(v);
		for(int t=0; t<d.cardT; ++t){
			SCIP_VAR * var;
			x_it[i].push_back(var);
			SCIPcreateVarBasic(scip, &(x_it[i][t]), ("x"+to_string(i)+to_string(t)).c_str(), 0, 1, 0, SCIP_VARTYPE_CONTINUOUS);
			SCIPaddVar(scip,x_it[i][t]);		
		}	
	}

	//contraintes sur x et y
	vector<vector<SCIP_CONS *> > cons_xy;
	for(int i=0; i<cardJ; ++i){
		vector<SCIP_CONS *> c;
		cons_xy.push_back(c);
		for(int t=0; t<cardT; ++t){
			SCIP_CONS * cons;
			cons_xy[i].push_back(cons);
			SCIPcreateConsLinear(scip, &cons_xy[i][t], ("cons_xy"+to_string(i)+to_string(t)).c_str(), 0, 0, 0, 0, 0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE);  
			SCIPaddCoefLinear(scip, cons_ct[i], x_it[i][t],1);		
			for(int l=0; l<L; ++l){
			//SCIPaddCoefLinear(scip, cons_ct[i], xt[i][j],-pente[i][j]);
			SCIPaddCoefLinear(scip, cons_ct[i], xtij[i][j],-pente[i][j]);
			SCIPaddCoefLinear(scip, cons_ct[i], pwd[i][j],-valbpt[i][j]+pente[i][j]*bpt[i][j]);
		}
		SCIPaddCons(scip, cons_ct[i]);
	}
	

}
