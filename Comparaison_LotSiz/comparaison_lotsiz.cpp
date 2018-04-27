#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <deque>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "data_struct.h"
//#include "Ordo_compact/ordo_comp.h"
#include "Lot_Sizing/Lot_Sizing_Dynamique/lotsizdyn2.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizingcom.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizcontcom.h"
#include "Lecteur_Fichiers/lecteur_pwd.h"
#include "Lecteur_Fichiers/lecteur_demande.h"
#include "Lecteur_Fichiers/generateur_pwd.h"
#include "Lecteur_Fichiers/generateur_demande.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizcontCPX.h"
#include "Lot_Sizing/Lot_Sizing_Compact/lotsizentCPX.h"
#include "Lecteur_Fichiers/fichier_sol.h"
#include "Lecteur_Fichiers/lecteur_sol.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
int infini = numeric_limits<int>::max();


int main(int argc, char *argv[]){
    /*if(argc != 6){ 
		cout << " Paramètres demandés : nombre périodes, nombre breakpoints, bornesup, stock maximal, numero instance" <<endl;
		return 0;
	}
    string cardT = argv[1];
	string nb_bp = argv[2];
    string bornesup = argv[3];
	string stockmax = argv[4];
    string num_instance = argv[5];*/


	vector<int> NBBP = {4,10};
    vector<int> S = {100,1000};
    vector<int> B = {100,2000};
    int cardT = 500;

    
    //int bornesup = 100;
    /*float tpsCPXcont = 0.0;
    float tpsCPXent = 0.0;
    float tpslotdyn = 0.0;
    float ecart = 0.0;*/
for(const auto& bornesup : B){
for(const auto& s : S){
for(const auto& nbbp : NBBP){

for(int i=1; i<11; ++i){
    //generateur_pwd(stoi(cardT), stoi(nb_bp), 5, stoi(bornesup));
    //generateur_demande(stoi(cardT),stoi(bornesup));

    string fichier_demande = "../Donnees/demande_"+to_string(cardT)+"_"+to_string(bornesup)+"_"+to_string(i);
    string fichier_pwd = "../Donnees/pwd_"+to_string(cardT)+"_"+to_string(nbbp)+"_"+to_string(bornesup)+"_"+to_string(i);

    //string fichier_demande = "../Donnees/demande_"+cardT+"_"+bornesup+"_"+num_instance;
    //string fichier_pwd = "../Donnees/pwd_"+cardT+"_"+nb_bp+"_"+bornesup+"_"+num_instance;


    data d;
    lecteur_demande(fichier_demande,d);
    //lecteur_demande("demande_20_100",d);
    lecteur_pwd(fichier_pwd,d);
    //lecteur_pwd("pwd_20_4_100",d);
    d.Q = s;
    d.s0 = 0;

    

// cout << d.bpt[0][1] << " "<< d.valbpt[0][1] << endl;
    /*vector<float> varia;
    auto start_time = chrono::steady_clock::now();
    lotsizcontcom(d,varia);
    auto end_time = chrono::steady_clock::now();
    cout << "temps contcom : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
    */
    float binfEnt = 0.0;
    float binfCont = 0.0;
    string statusEnt, statusCont;
    float tpsCPXent, tpsCPXcont;

    auto start_time = chrono::steady_clock::now();
    float solDyn = lotsizdyn(d, 1);
    auto end_time = chrono::steady_clock::now();
    float tpslotdyn = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
    //cout << "temps dyn : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
    

   /* start_time = chrono::steady_clock::now();
    lotsizcomp(d, 1);
    end_time = chrono::steady_clock::now();
    cout << "temps com : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
*/
    //start_time = chrono::steady_clock::now();
    float solent = lotsizentCPX(d, 1,binfEnt,statusEnt, tpsCPXent);
    //end_time = chrono::steady_clock::now();
    //float tpsCPXent = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
    //cout << "temps com : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
    
    //start_time = chrono::steady_clock::now();
    float solcont = lotsizcontCPX(d,binfCont,statusCont,tpsCPXcont);
    //end_time = chrono::steady_clock::now();
    //float tpsCPXcont = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
    
    //float ecart = (solent - solcont)*100/solcont;
//}
/*
    cout << "temps moyen lotsizdyn : " << tpslotdyn/10 << endl;
    cout << "temps moyen CPLEX ent : " << tpsCPXent/10 << endl;
    cout << "temps moyen CPLEX cont : " << tpsCPXcont/10 << endl;
    cout << "ecart moyen : " << ecart/10 << endl;*/

    fichier_sol(cardT,nbbp,bornesup,d.Q,i,tpslotdyn,solDyn,tpsCPXent,solent,tpsCPXcont,solcont,binfEnt,binfCont,statusEnt,statusCont);
    //fichier_sol(d.cardT,nbbp,bornesup,d.Q,i,tpslotdyn,solDyn,tpsCPXent,solent,tpsCPXcont,solcont);
}
}
}
}
    //lecteur_sol("../Log_sol/"+cardT+"_"+nb_bp+"_"+bornesup+"_"+stockmax);
    return 0;
}