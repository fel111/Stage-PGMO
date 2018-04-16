#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <algorithm>
//#include <chrono>
//#include <deque>
//#include "scip/scip.h"
//#include "scip/scipdefplugins.h"
#include "data_struct.h"
#include "Ordo_compact/ordo_comp.h"
#include "Ordo_compact/ordo_cplex.h"
#include "Lot_Sizing/Lot_Sizing_Dynamique/lotsizdyn2.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizentCPX.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizcontcom.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizcontCPX.h"
#include "Modele_compact/compact.h"
#include "Modele_compact/modele_entier_cplex.h"
#include "Lecteur_Fichiers/lecteur_taches.h"
#include "Lecteur_Fichiers/lecteur_pwd.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
int infini = numeric_limits<int>::max();


int main(){
	data d;
	// init test

	d.cardJ = 2;
	d.cardM = 1;
	d.cardR = 2;
	d.s0 = 0;
	d.cardT = 5; 
	d.Q = 100;
	for(int i=0;i<d.cardT;++i){
		vector<float> temp_pente = {2.0,1.0};
		vector<int> temp_bpt = {0,20,infini};
		vector<float> temp_valbpt = {0.0,40.0};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(2);
	}

	d.Ckr = {{100000.0,50.0}};
	d.Dk = {0.0,0.0}; // ne converge pas avec != 0

	lecteur_taches("../Donnees/taches_5",d);
	//lecteur_pwd("../Donnees/pwd_20_4_2000_1",d);
	data dinit = d;
	float solOrdo, solLotSiz;
	vector<float> rt;

	cout << "MODELE COMPACT CPLEX :"<<endl;
	modele_entier_cplex(d);

	cout << "DECOMPOSITION"<< endl;
	solOrdo = ordo_cplex(d);
	//solOrdo = ordo(d);
	cout << "sol ordocplex : "<< solOrdo << endl;
	
	solLotSiz = lotsizcontCPX(d, rt);
	//solLotSiz = lotsizcontcom(d, rt);
	cout << "sol lotsizcplex : "<< solLotSiz << endl << endl;

	/*for(int t=0; t<d.cardT; ++t){
		cout << "bpt et valbpt [0] avant modif: " << d.bpt[t][0] << " " << d.valbpt[t][0] << endl;
	}*/
	modifPWL(d, rt);
	/*for(int t=0; t<d.cardT; ++t){
		cout << "bpt et valbpt [0] après modif: " << d.bpt[t][0] << " " << d.valbpt[t][0] << endl;
	}*/


	int i=0;



	while(i<1){
	//while(solOrdo != solLotSiz){


		solOrdo = ordo_cplex(d);
		//solOrdo = ordo(d);
		cout << "sol ordocplex : "<< solOrdo << endl;
		//cout << "ORDO REALISE" << endl;
		//for(int i=0; i<d.cardT; ++i) cout << d.dt[i] <<" ";
		//cout << endl;
		
		initPWD(dinit,d);
		//cout << "LOTSIZING CONT COMP"<<endl;
		solLotSiz = lotsizcontCPX(d, rt);
		//solLotSiz = lotsizcontcom(d, rt);
		cout << "sol lotsizcplex : "<< solLotSiz << endl;

		//cout << "LOTSIZING DYN"<<endl;
		//vector<int> rt = lotsizdyn(d,1);
		
		cout << endl;
		//lotsizdyn(d,2);

		
		//cout << "MODIFICATION PWL" << endl;
		modifPWL(d, rt);
		
		++i;
	}




	

	//cout << "LOTSIZING COMP"<<endl;
	//lotsizcomp(d);


	return 0;
}








