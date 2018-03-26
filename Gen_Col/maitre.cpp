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
	int cardT, cardJ, cardM, cardR;
	vector<int> pj;
	vector<vector<float> > cjr;
	vector<vector<float> > Ckr;
	vector<float> Dk;
	vector<vector<float> > Djk;
	vector<vector<int> > ri;
	vector<vector<int> > di;
};



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


vector<vector<float> > firstSol(data d){
	<vector<float> cl (d.cardT, 0.0);
	vector<vector<float> > ress = d.Ckr;
	vector<vector<vector<int> > > yjkt (d.cardJ, vector<vector<int> > (d.cardM, vector<int> (d.cardT, 0)));
	for(int j=0; j<d.cardJ; ++j){
		//int fixed = false;	
		for(int t=0; t<d.cardT; ++t){
			for(int k=0; k<d.cardM; ++k){
				bool test = true;	
				for(int r=0; r<d.cardR; ++r){
					if((d.cjr[j][r]<=ress[k][r])&&(d.ri[j]>=t)&&(t+d.pj[j]<=d.dj[j])&&(fixed==false)){
						ress[k][r] = ress[k][r]-d.cjr[j][r];
						fixed = true;
						cl[fixed] += p_t(t,d.Djk[j][k],d);
						yjkt[j][k][t] = 1;	
					}
				}
			}
			
		}
	}
	
}


int main(){
	
	
	
	// MODELE



	SCIP * scip;
	SCIPcreate(&scip);
	SCIPincludeDefaultPlugins(scip);
	SCIPcreateProb(scip, "PMR", NULL, NULL, NULL, NULL, NULL, NULL, NULL);

	

	

}
