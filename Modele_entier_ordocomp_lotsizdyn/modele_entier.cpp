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
	if(argc != 3){
		cout << "erreur param : fichier instance, fichier parametre" << endl;
		return 0;
	}
	
	string instance = argv[1];
	string param = argv[2];

	d.cardR = 2;
	d.s0 = 0;
	d.Q = 100;
	d.cardM = 1;


	//lecteur_taches("../Donnees/taches_"+to_string(d.cardT),d, true);
	lecteur_taches_EnergSchedInst("../Donnees/dataSchedulingInstances/inst30-1_25-0_1-0_8-10.dat",d);

	//cout << "cardT : " << d.cardT << endl;

	for(int i=0;i<d.cardT;++i){
		vector<float> temp_pente = {1.0,0.5,2.0,1.5};
		vector<float> temp_bpt = {0.0,4.0,8.0,100.0,infini};
		vector<float> temp_valbpt = {0.0,4.0,6.0,190.0};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(4);
	}
	vector<vector<float> > ress;
	for(int i=0;i<d.cardM;++i){
		//vector<float> r = {(d.cardJ*4.0/d.cardM),(d.cardJ*4000.0/d.cardM)}; // s'assure qu'il y ai assez de ram/cpu
		vector<float> r = {d.cardJ, d.cardJ};
		//cout << r[0] << " " << r[1] << endl;
		ress.push_back(r);
	}
	d.Ckr = ress;
	vector<float> Dk (d.cardM, 0.0); // ne converge pas avec != 0
	d.Dk = Dk; 

	//lecteur_taches("../Donnees/taches_"+to_string(d.cardT),d);
	//lecteur_pwd("../Donnees/pwd_20_4_2000_4",d);
	lecteur_param("../Param/param1.txt",p);
	data dinit = d;
	float solOrdo, solLotSizDyn, solLotSizMip, solLotSizPLNE, solComp, solMip, solPLNE, solDyn, tpsBoucleMip, tpsBouclePLNE, tpsBoucleDyn, tpsComp, tpsRelax, tpsLS, borneinfComp,borneinfRelax,borneinftemp;
	string statusComp, statusRelax, statusOrdo, statusLS;
	vector<float> rt;
	int boucle;
	vector<float> solrelax;
	float tpsOrdo;
	float tpsTotOrdo = 0.0;
	float tpsTotLS = 0.0;




	// MODELE COMPACT CPLEX
	//auto start_time = chrono::steady_clock::now();
	solComp = modele_entier_cplex(d,p,tpsComp,borneinf,statusComp);
	//auto end_time = chrono::steady_clock::now();
	//tpsComp = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	cout << "COMPACT		: tps sol    " << tpsComp << " " << solComp << " " << statusComp << " " << borneinfComp << endl;

	
	if(p.boucle_relaxation==1){
		//start_time = chrono::steady_clock::now();
		relaxation_modele_entier_cplex(d,solrelax,p,tpsRelax, borneinfRelax,statusRelax);
		//end_time = chrono::steady_clock::now();
		//tpsRelax = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
		cout << "RELAX tps= "<<tpsRelax<<" status="<<statusRelax<<" borneinf="<<borneinfRelax<<endl;
	}




	// BOUCLE ORDO + LOTSIZING PL (PWD)
	auto start_time = chrono::steady_clock::now();
	if(p.boucle_relaxation==1){
		modifPWL(d, solrelax);
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		initPWD(dinit,d);
		boucle = p.nb_iter_max_boucle - 1;
	} 
	else{
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		boucle = p.nb_iter_max_boucle ;
	}
	tpsTotOrdo += tpsOrdo;
	solLotSizMip = lotsizcontCPX(d, p, rt, tpsLS,borneinftemp,statusLS);
	tpsTotLS += tpsLS;
	boucle = p.nb_iter_max_boucle - 1;
	while((solOrdo != solLotSizMip)&&(boucle>0)){
		modifPWL(d, rt);
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		tpsTotOrdo += tpsOrdo;
		//cout << "solOrdo : "<<solOrdo<<endl;
		initPWD(dinit,d);
		solLotSizMip = lotsizcontCPX(d, p, rt, tpsLS,borneinftemp,statusLS);
		tpsTotLS += tpsLS;
		//cout << "solLS : "<<solLotSizMip<<endl;
		--boucle;
	}
	end_time = chrono::steady_clock::now();
	tpsBoucleMip = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solMip = solLotSizMip;

	cout << "BOUCLE PL (PWD)  tpsBoucle=" << tpsBoucleMip << " tpsOrdo="<<tpsTotOrdo<<" tpsLS="<<tpsTotLS<<" bornesup=" << solMip;
	if (boucle==0) cout << " max iter atteint " << endl;
	else cout << " ordo = lotsizing" << endl;
	d = dinit;




	// BOUCLE ORDO + LOTSIZING PLNE
	start_time = chrono::steady_clock::now();
	if(p.boucle_relaxation==1){
		modifPWL(d, solrelax);
		solOrdo = ordo_cplex(d,p);
		initPWD(dinit,d);
	} 
	else solOrdo = ordo_cplex(d,p);
	solLotSizPLNE = lotsizentCPX(d,1,rt, p);
	boucle = p.nb_iter_max_boucle - 1;
	while((solOrdo != solLotSizPLNE)&&(boucle>0)){
	//while((!conv)&&(boucle>0)){
		modifPWL(d, rt);
		solOrdo = ordo_cplex(d,p);
		//cout << "sol plne ordo : " << solOrdo << endl;
		//if(solOrdo == solOrdoTemp) conv = true;
		//else{
		initPWD(dinit,d);
		solLotSizPLNE = lotsizentCPX(d,1,rt, p);
		//cout << "sol plne ls : " << solLotSizPLNE << endl;
		
		//solOrdoTemp = solOrdo;
		//}
		--boucle;
	}
	end_time = chrono::steady_clock::now();
	tpsBouclePLNE = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solPLNE = solOrdo;

	cout << "BOUCLE PLNE	: tps sol    " << tpsBouclePLNE << " " << solPLNE;
	if (boucle==0) cout << " max iter atteint " << endl;
	else cout << " ordo = lotsizing" << endl;

	d = dinit;




	// BOUCLE ORDO + LOTSIZING DYNAMIQUE
	start_time = chrono::steady_clock::now();
	if(p.boucle_relaxation==1){
		modifPWL(d, solrelax);
		solOrdo = ordo_cplex(d,p);
		initPWD(dinit,d);
	} 
	else solOrdo = ordo_cplex(d,p);
	solLotSizDyn = lotsizdyn(d,1,rt);
	boucle = p.nb_iter_max_boucle - 1;

	//while((!conv)&&(boucle>0)){
	while((solOrdo != solLotSizDyn)&&(boucle>0)){
		modifPWL(d, rt);
		solOrdo = ordo_cplex(d,p);
		//cout <<  "sol ordo lsdyn : " << solOrdo <<endl;
		//if(solOrdo == solOrdoTemp) conv = true;
		//else{
		initPWD(dinit,d);
		solLotSizDyn = lotsizdyn(d,1,rt);
			//cout << "sol lsdyn : " << solLotSizDyn << endl;
			
			//solOrdoTemp = solOrdo;
		//}
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








