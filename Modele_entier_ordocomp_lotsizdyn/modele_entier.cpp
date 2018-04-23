#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <algorithm>
#include <chrono>
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
#include "Lecteur_Fichiers/lecteur_param.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
float infini = numeric_limits<float>::max();


int main(int argc, char* argv[]){
	data d;
	param p;
	// init test
	if(argc != 4){
		cout << "erreur param : nombre périodes, nombre machines, nombre taches" << endl;
		return 0;
	}
	
	d.cardT = stoi(argv[1]);
	d.cardM = stoi(argv[2]);
	d.cardJ = stoi(argv[3]);
	
	d.cardR = 2;
	d.s0 = 0;
	d.Q = 100;

	for(int i=0;i<d.cardT;++i){
		vector<float> temp_pente = {1.0,1.2,4.8,3.1,1.6,4.7,1.7,2.0,3.8,2.7};
		vector<float> temp_bpt = {0.0,200.0,400.0,600.0,800.0,1000.0,1200.0,1400.0,1600.0,1800.0,infini};
		vector<float> temp_valbpt = {0.0,200.0,440.0,1400.0,2020.0,2340.0,3280.0,3620.0,4020.0,4780.0};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(10);
	}
	vector<vector<float> > ress;
	for(int i=0;i<d.cardM;++i){
		vector<float> r = {(d.cardJ*4.0/d.cardM),(d.cardJ*4000.0/d.cardM)}; // s'assure qu'il y ai assez de ram/cpu
		//cout << r[0] << " " << r[1] << endl;
		ress.push_back(r);
	}
	d.Ckr = ress;
	vector<float> Dk (d.cardM, 0.0); // ne converge pas avec != 0
	d.Dk = Dk; 

	lecteur_taches("../Donnees/taches_20",d,true);
	//lecteur_pwd("../Donnees/pwd_20_4_2000_4",d);
	lecteur_param("../Param/param1.txt",p);
	data dinit = d;
	float solOrdo, solLotSizDyn, solLotSizMip, solLotSizPLNE, solComp, solMip, solPLNE, solDyn, tpsBoucleMip, tpsBouclePLNE, tpsBoucleDyn, tpsComp;
	vector<float> rt;
	int boucle;



	// MODELE COMPACT CPLEX
	auto start_time = chrono::steady_clock::now();
	solComp = modele_entier_cplex(d,p);
	auto end_time = chrono::steady_clock::now();
	tpsComp = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	cout << "COMPACT     : tps sol    " << tpsComp << " " << solComp << endl;

	// BOUCLE ORDO + LOTSIZING PL (PWD)
	//cout << "relax : " << p.boucle_relaxation << endl;
	start_time = chrono::steady_clock::now();
	if(p.boucle_relaxation==1){
		float solrelax = relaxation_modele_entier_cplex(d,rt,p);
		//cout << "solrelax : " << solrelax << endl;
		modifPWL(d, rt);
		solOrdo = ordo_cplex(d,p);
		initPWD(dinit,d);
	} 
	else solOrdo = ordo_cplex(d,p);
	solLotSizMip = lotsizcontCPX(d, rt, p);
	boucle = 2;
	while((solOrdo != solLotSizMip)&&(boucle>0)){
		modifPWL(d, rt);
		//cout << "------------------ORDO ----------------------"<<endl;
		solOrdo = ordo_cplex(d,p);
		//cout << "solOrdo : "<<solOrdo<<endl;
		initPWD(dinit,d);
		//cout << "------------------LOTSIZING ----------------------"<<endl;
		solLotSizMip = lotsizcontCPX(d, rt, p);
		//cout << "solLS : "<<solLotSizMip<<endl;
		--boucle;
		
	}
	end_time = chrono::steady_clock::now();
	tpsBoucleMip = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solMip = solLotSizMip;

	cout << "BOUCLE PL (PWD)  : tps sol    " << tpsBoucleMip << " " << solMip;
	if (boucle==0) cout << " max iter atteint " << endl;
	else cout << " ordo = lotsizing" << endl;
	d = dinit;

	// BOUCLE ORDO + LOTSIZING PLNE
	start_time = chrono::steady_clock::now();
	float solOrdoTemp = ordo_cplex(d,p);
	bool conv = false;
	solLotSizPLNE = lotsizentCPX(d,1,rt, p);
	modifPWL(d, rt);
	boucle = 2;


	while((!conv)&&(boucle>0)){
		solOrdo = ordo_cplex(d,p);
		if(solOrdo == solOrdoTemp) conv = true;
		else{
			initPWD(dinit,d);
			solLotSizPLNE = lotsizentCPX(d,1,rt, p);
			modifPWL(d, rt);
			solOrdoTemp = solOrdo;
		}
		--boucle;
	}
	end_time = chrono::steady_clock::now();
	tpsBouclePLNE = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solPLNE = solOrdo;

	cout << "BOUCLE PLNE : tps sol    " << tpsBouclePLNE << " " << solPLNE;
	if (boucle==0) cout << " max iter atteint " << endl;
	else cout << " ordo = lotsizing" << endl;

	d = dinit;

	// BOUCLE ORDO + LOTSIZING DYNAMIQUE
	start_time = chrono::steady_clock::now();
	solOrdoTemp = ordo_cplex(d,p);
	conv = false;
	solLotSizDyn = lotsizdyn(d,1,rt);
	modifPWL(d, rt);
	boucle = 2;


	while((!conv)&&(boucle>0)){
		solOrdo = ordo_cplex(d,p);
		if(solOrdo == solOrdoTemp) conv = true;
		else{
			initPWD(dinit,d);
			solLotSizDyn = lotsizdyn(d,1,rt);
			modifPWL(d, rt);
			solOrdoTemp = solOrdo;
		}
		--boucle;
	}
	end_time = chrono::steady_clock::now();
	tpsBoucleDyn = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solDyn = solOrdo;

	cout << "BOUCLE DYN  : tps sol    " << tpsBoucleDyn << " " << solDyn;
	if (boucle==0) cout << " max iter atteint " << endl;
	else cout << " ordo = lotsizing" << endl;
	return 0;
}








