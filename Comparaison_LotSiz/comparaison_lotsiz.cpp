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
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
int infini = numeric_limits<int>::max();


int main(int argc, char *argv[]){
    if(argc != 5){ 
		cout << " Paramètres demandés : nombre périodes, nombre breakpoints, bornesup, stock maximal" <<endl;
		return 0;
	}
    string cardT = argv[1];
	string nb_bp = argv[2];
    string bornesup = argv[3];
	string stockmax = argv[4];
	
    float tpsCPXcont = 0.0;
    float tpsCPXent = 0.0;
    float tpslotdyn = 0.0;
    float ecart = 0.0;

    for(int i=0; i<10; ++i){

        generateur_pwd(stoi(cardT), stoi(nb_bp), 5, stoi(bornesup));
        generateur_demande(stoi(cardT),stoi(bornesup));

        string fichier_demande = "demande_"+cardT+"_"+bornesup;
        string fichier_pwd = "pwd_"+cardT+"_"+nb_bp+"_"+bornesup;



        data d;
        lecteur_demande(fichier_demande,d);
        //lecteur_demande("demande_20_100",d);
        lecteur_pwd(fichier_pwd,d);
        //lecteur_pwd("pwd_20_4_100",d);
        d.Q = stoi(stockmax);
        d.s0 = 0;

        

        /*for(int i=0; i<d.cardT; ++i){
            cout << d.dt[i] << endl;
        }*/

    // cout << d.bpt[0][1] << " "<< d.valbpt[0][1] << endl;
        /*vector<float> varia;
        auto start_time = chrono::high_resolution_clock::now();
        lotsizcontcom(d,varia);
        auto end_time = chrono::high_resolution_clock::now();
        cout << "temps contcom : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
        */

        auto start_time = chrono::high_resolution_clock::now();
        lotsizdyn(d, 1);
        auto end_time = chrono::high_resolution_clock::now();
        tpslotdyn += chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
        //cout << "temps dyn : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
        

        /*start_time = chrono::high_resolution_clock::now();
        lotsizcomp(d, 1);
        end_time = chrono::high_resolution_clock::now();
        cout << "temps com : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
    */
        start_time = chrono::high_resolution_clock::now();
        float solent = lotsizentcplex(d, 1);
        end_time = chrono::high_resolution_clock::now();
        tpsCPXent += chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
        //cout << "temps com : " << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0 << endl<<endl;
        start_time = chrono::high_resolution_clock::now();
        float solcont = lotsizcontcplex(d);
        end_time = chrono::high_resolution_clock::now();
        tpsCPXcont += chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
        
        ecart += (solent - solcont)*100/solcont;
    }

    cout << "temps moyen lotsizdyn : " << tpslotdyn/10 << endl;
    cout << "temps moyen CPLEX ent : " << tpsCPXent/10 << endl;
    cout << "temps moyen CPLEX cont : " << tpsCPXcont/10 << endl;
    cout << "ecart moyen : " << ecart/10 << endl;

    return 0;
}