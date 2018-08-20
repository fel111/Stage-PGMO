#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "../struct.h"
#include "lecteur_param.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;


int lecteur_param(string file, param &p, data &d){
	
    ifstream fichier(file, ios::in); //ouverture du fichier
    if(fichier){  // si l'ouverture fonctionne
        string ligne;
        string par;
        string t;

        while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
            istringstream iss (ligne);
            iss >> par >> t;
            if (par == "taches_avec_fenetre_temps") { p.taches_avec_fenetre_temps = stoi(t); }
            //else if (par == "boucle_relaxation") { p.boucle_relaxation = stoi(t); }
            else if (par == "aff_log_ordo_cplex") { p.aff_log_ordo_cplex = stoi(t); } 
            //else if (par == "aff_log_lotsizingcont_cplex") { p.aff_log_lotsizingcont_cplex = stoi(t); }
            //else if (par == "aff_log_lotsizingent_cplex") { p.aff_log_lotsizingent_cplex = stoi(t); }
            else if (par == "aff_log_compact_cplex") { p.aff_log_compact_cplex = stoi(t); }
            else if (par == "nb_threads_cplex") { p.nb_threads_cplex = stoi(t); }
            else if (par == "time_limit_ordo") { p.time_limit_ordo = stoi(t); }
            //else if (par == "time_limit_lotsizingcont") { p.time_limit_lotsizingcont = stoi(t); }
            //else if (par == "time_limit_lotsizingent") { p.time_limit_lotsizingent = stoi(t); }
            //else if (par == "time_limit_lotsizingdyn") { p.time_limit_lotsizingdyn = stoi(t); }
            else if (par == "time_limit_compact") { p.time_limit_compact = stoi(t); }
            //else if (par == "time_limit_relax") { p.time_limit_relax = stoi(t); }
            //else if (par == "nb_iter_max_boucle") { p.nb_iter_max_boucle = stoi(t); }
            //else if (par == "choix_dt_ls") { p.choix_dt_ls = stoi(t); }
            else if (par == "type_cap") { p.type_cap = stoi(t); }
            else if (par == "val_cap") { p.val_cap = stof(t); }
            else if (par == "time_limit_gencol") { p.time_limit_gencol = stoi(t); }
            else if (par == "ensemble_multiple") p.ensemble_multiple = stoi(t);
            else if (par == "ajout_breakpoint_voisin") p.ajout_breakpoint_voisin = stoi(t);
        }
        if((p.type_cap == 1)&&(p.val_cap <1)){
            cout << "type capacite : abs  mais capacite < 1   erreur"<<endl;
            return 0;
        }
	    if((p.type_cap == 2)&&((p.val_cap <0)||(p.val_cap>1))){
            cout << "type capacite : rel  mais capacite != {0,1}   erreur"<<endl;
            return 0;
        }
	

        return 1;
    }
    else{
        cout << "erreur lecture fichier param" << endl;
        return 0;
    }
}
