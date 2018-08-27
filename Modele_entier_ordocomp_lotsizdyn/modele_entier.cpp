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
#include "../struct.h"
//#include "Ordo_compact/ordo_comp.h"
#include "../Ordo_compact/ordo_cplex.h"
#include "../Lot_Sizing/Lot_Sizing_Dynamique/lotsizdyn2.h"
//#include "Lot_Sizing/Lot_Sizing_Compact/lotsizentCPX.h"
//#include "Lot_Sizing/Lot_Sizing_Compact/lotsizcontcom.h"
//#include "Lot_Sizing/Lot_Sizing_Compact/lotsizcontCPX.h"
//#include "Modele_compact/compact.h"
#include "../Modele_compact/modele_entier_cplex.h"
#include "../Lecteur_Fichiers/lecteur_taches.h"
#include "../Lecteur_Fichiers/lecteur_pwl.h"
#include "../Lecteur_Fichiers/lecteur_param.h"
#include "../Gen_Col_noStock/genColNoStock.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
//float infini = numeric_limits<float>::max();


int main(int argc, char* argv[]){
	data d;
	param p;
	// init test
	if(argc != 4){
		cout << "erreur param : fichier instance_taches, fichier instance_pwl, fichier parametre" << endl;
		return 0;
	}
	
	string instance_tache = argv[1];
	string instance_pwl = argv[2];
	string parametre = argv[3];

	/*d.cardR = 2;
	d.s0 = 0;
	//d.Q = 100;
	d.cardM = 1;*/

	if(lecteur_param("Param/"+parametre,p,d) == 0) return 0;
	
	lecteur_taches_EnergSchedInst("Donnees/dataSchedulingInstances_newname/"+instance_tache,d);

	lecteur_pwl("Donnees/pwl/"+instance_pwl,d);

	initCap(d,p);
	//lecteur_taches_EnergSchedInst("../Donnees/dataSchedulingInstances/inst60-1_25-0_1-3_2-1.dat",d);

	//cout << "cardT : " << d.cardT << endl;

	/*for(int i=0;i<d.cardT;++i){
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
	d.Dk = Dk; */

	//lecteur_taches("../Donnees/taches_"+to_string(d.cardT),d);
	//lecteur_pwd("../Donnees/pwd_20_4_2000_4",d);
	
	data dinit = d;
	float solOrdo, solLS, solComp, solMip, solPLNE, solDyn, tpsBoucleMip, tpsBouclePLNE, tpsBoucleDyn, tpsComp, tpsRelax, tpsLS, borneinfComp,borneinfRelax,borneinftemp;
	string statusComp, statusRelax, statusOrdo, statusLS;
	vector<float> rt, demande;
	int boucle;
	vector<float> solrelax;
	float tpsOrdo;
	float tpsTotOrdo = 0.0;
	float tpsTotLS = 0.0;
	//float timeBoucle = p.time_limit_compact / p.nb_iter_max_boucle;
	float timeRab = 0.0;
	float timeRabRelax = 0.0;
	int cptBoucle = 0;

	//cout << "instance="<<instance<<" parametre="<<parametre<<" nbPeriodes="<<d.cardT<<" nbTaches="<<d.cardJ<<" nbMachines="<<d.cardM<<" borneBatterie="<<d.Q<<" ratioDemandeTot/Capacite="<<d.Q/d.consoTot<<" limiteNbBoucles=-1 limiteTpsOrdo=300 limiteThreads="<<p.nb_threads_cplex<<endl;
	// MODELE COMPACT CPLEX
	//auto start_time = chrono::steady_clock::now();
	//solComp = modele_entier_cplex(d,p,tpsComp,borneinfComp,statusComp);
	//auto end_time = chrono::steady_clock::now();
	//tpsComp = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	//cout << "COMPACT tps=" << tpsComp << " statut=" << statusComp << " bornesup=" << solComp  << " borneinf=" << borneinfComp << endl;

	
	/*if(p.boucle_relaxation==1){
		p.time_limit_relax = timeBoucle;
		//start_time = chrono::steady_clock::now();
		float borneSupRelax = relaxation_modele_entier_cplex(d,solrelax,p,tpsRelax, borneinfRelax,statusRelax);
		//end_time = chrono::steady_clock::now();
		//tpsRelax = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
		cout << "INIT=relax tps= "<<tpsRelax<<" statut="<<statusRelax<<" bornesup="<<borneSupRelax<<" borneinf="<<borneinfRelax<<endl;;
		timeRabRelax = timeBoucle-tpsRelax;
	}
	else{
		cout << "INIT=none"<<endl;
	}*/

	/*for(int i=0; i<d.cardT;++i){
		cout << "relax["<<i<<"] = "<<solrelax[i] <<endl;
	}*/


	// BOUCLE ORDO + LOTSIZING PL (PWD)
	/*auto start_time = chrono::steady_clock::now();
	if(p.boucle_relaxation==1){
		p.time_limit_ordo = timeBoucle/2+timeRabRelax;
		modifPWL(d, solrelax);
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		
		initPWD(dinit,d);
		boucle = p.nb_iter_max_boucle - 1;
	} 
	else{
		p.time_limit_ordo = timeBoucle/2;
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		boucle = p.nb_iter_max_boucle ;
	}
	//cout << "solOrdo : "<<solOrdo<< endl;
	tpsTotOrdo += tpsOrdo;
	timeRab = p.time_limit_ordo - tpsOrdo;
	p.time_limit_lotsizingcont = timeBoucle/2 + timeRab;
	solLS = lotsizcontCPX(d, p, rt, tpsLS,borneinftemp,statusLS);
	//cout << "solLS : "<<solLS<< endl;
	timeRab = p.time_limit_lotsizingcont - tpsLS;
	tpsTotLS += tpsLS;
	++cptBoucle;
	--boucle;
	while((solOrdo != solLS)&&(boucle>0)){
		modifPWL(d, rt);
		p.time_limit_ordo = timeBoucle/2 + timeRab;
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		//cout << "solOrdo : "<<solOrdo<< endl;
		timeRab = p.time_limit_ordo - tpsOrdo;
		tpsTotOrdo += tpsOrdo;
		//cout << "solOrdo : "<<solOrdo<<endl;
		initPWD(dinit,d);
		p.time_limit_lotsizingcont = timeBoucle/2 + timeRab;
		solLS = lotsizcontCPX(d, p, rt, tpsLS,borneinftemp,statusLS);
		//cout << "solLS : "<<solLS<< endl;
		timeRab = p.time_limit_lotsizingcont - tpsLS;
		tpsTotLS += tpsLS;
		//cout << "solLS : "<<solLotSizMip<<endl;
		--boucle;
		++cptBoucle;
	}
	auto end_time = chrono::steady_clock::now();
	tpsBoucleMip = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solMip = solLS;

	cout << "BOUCLE PL (PWD)  tpsBoucle=" << tpsBoucleMip << " tpsOrdo="<<tpsTotOrdo<<" tpsLS="<<tpsTotLS<<" bornesup=" << solMip <<" nbBoucle="<<cptBoucle;
	if (solOrdo==solLS) cout << " arret=CONV" << endl;
	else cout << " arret=MAXITER" << endl;
	d = dinit;
	tpsTotOrdo = 0.0;
	tpsTotLS = 0.0;
	cptBoucle = 0;*/
	float _delete;
	bool solOrdoOk = false;
	// BOUCLE ORDO_GENCOL + LOTSIZING DYNAMIQUE
	auto start_time = chrono::steady_clock::now();
	solOrdo = genColNoStock(d,p,tpsOrdo,demande,_delete,statusOrdo);
	tpsTotOrdo += tpsOrdo;
	cout << "solordocplex : " << ordo_cplex(d,p,tpsBoucleMip,demande,statusComp,tpsBoucleDyn)<<endl;
	cout << "solordo : " << solOrdo << endl;
	if(solOrdo != -1){
		solOrdoOk = true;
		solLS = lotsizdyn(d,p,rt,tpsLS,demande);
		cout << "solLS : " << solLS << endl;
		tpsTotLS += tpsLS;
		++cptBoucle;
	}
	while((solOrdo != solLS) && (tpsTotOrdo < 7200) && (solOrdoOk)){
		modifPWL(d, rt);
		solOrdo = genColNoStock(d,p,tpsOrdo,demande,_delete,statusOrdo);
		cout << "solordocplex : " << ordo_cplex(d,p,tpsBoucleMip,demande,statusComp,tpsBoucleDyn)<<endl;
		cout << "solordo : " << solOrdo << endl;
		tpsTotOrdo += tpsOrdo;
		if(solOrdo == -1) solOrdoOk = false;
		else{
			initPWD(dinit,d);
			solLS = lotsizdyn(d,p,rt,tpsLS,demande);
			cout << "solLS : " << solLS << endl;
			tpsTotLS += tpsLS;
			++cptBoucle;
		}
	}
	auto end_time = chrono::steady_clock::now();
	tpsBoucleDyn = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	cout << "solOrdo/LS : " << solLS << endl;

	// BOUCLE ORDO + LOTSIZING DYNAMIQUE
	/*auto start_time = chrono::steady_clock::now();
	/*if(p.boucle_relaxation==1){
		modifPWL(d, solrelax);
		p.time_limit_ordo = timeBoucle/2+timeRabRelax;
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		initPWD(dinit,d);
		boucle = p.nb_iter_max_boucle - 1;
	} 
	else{* /
		p.time_limit_ordo = 300;
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		//boucle = p.nb_iter_max_boucle ;
	//}
	//cout << "solOrdo : "<<solOrdo<< endl;
	tpsTotOrdo += tpsOrdo;
	//timeRab = p.time_limit_ordo - tpsOrdo;
	p.time_limit_lotsizingdyn = 300;
	solLS = lotsizdyn(d,p,rt,tpsLS);
	//cout << "solLS : "<<solLS<< endl;
	//for(int i=0; i<rt.size(); ++i) cout << rt[i] << " ";
	//cout << endl;
	//timeRab = p.time_limit_lotsizingcont - tpsLS;
	tpsTotLS += tpsLS;
	++cptBoucle;
	//--boucle;
	//while((!conv)&&(boucle>0)){
	//while((solOrdo != solLS)&&(boucle>0)){
	while(solOrdo != solLS){
		modifPWL(d, rt);
		//p.time_limit_ordo = timeBoucle/2 + timeRab;
		solOrdo = ordo_cplex(d,p,tpsOrdo,statusOrdo);
		//cout << "solOrdo : "<<solOrdo<< endl;
		//timeRab = p.time_limit_ordo - tpsOrdo;
		tpsTotOrdo += tpsOrdo;
		//cout << "solOrdo : "<<solOrdo<<endl;
		initPWD(dinit,d);
		//p.time_limit_lotsizingdyn = timeBoucle/2 + timeRab;
		solLS = lotsizdyn(d,p,rt,tpsLS);
		//cout << "solLS : "<<solLS<< endl;
		//timeRab = p.time_limit_lotsizingcont - tpsLS;
		tpsTotLS += tpsLS;
			//cout << "sol lsdyn : " << solLotSizDyn << endl;
			
			//solOrdoTemp = solOrdo;
		//}
		//--boucle;
		++cptBoucle;
	}
	auto end_time = chrono::steady_clock::now();
	tpsBoucleDyn = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
	solDyn = solOrdo;
	cout << "BOUCLE Dyn  tpsBoucle=" << tpsBoucleDyn << " tpsOrdo="<<tpsTotOrdo<<" tpsLS="<<tpsTotLS<<" bornesup=" << solDyn << " nbBoucle="<<cptBoucle<<" statutDernierOrdo="<<statusOrdo<<endl;
	//if (solOrdo==solLS) cout << " arret=CONV" << endl;
	//else cout << " arret=MAXITER" << endl;

	//cout << endl << endl;
	*/
	return 0;
}








