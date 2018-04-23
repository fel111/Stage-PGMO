#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "data_struct.h"
#include "lecteur_param.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;


void lecteur_param(string file, param &p){
	
    ifstream fichier(file, ios::in); //ouverture du fichier
    if(fichier){  // si l'ouverture fonctionne
        string ligne;
        string par;
        int t;
        while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
            istringstream iss (ligne);
            iss >> par >> t;
            if (par == "taches_avec_fenetre_temps") { p.taches_avec_fenetre_temps = t; }
            else if (par == "boucle_relaxation") { p.boucle_relaxation = t; }
            else if (par == "aff_log_ordo_cplex") { p.aff_log_ordo_cplex = t; } 
            else if (par == "aff_log_lotsizingcont_cplex") { p.aff_log_lotsizingcont_cplex = t; }
            else if (par == "aff_log_lotsizingent_cplex") { p.aff_log_lotsizingent_cplex = t; }
            else if (par == "aff_log_compact_cplex") { p.aff_log_compact_cplex = t; }
            else if (par == "nb_threads_cplex") { p.nb_threads_cplex = t; }
            else if (par == "time_limit_ordo") { p.time_limit_ordo = t; }
            else if (par == "time_limit_lotsizingcont") { p.time_limit_lotsizingcont = t; }
            else if (par == "time_limit_lotsizingent") { p.time_limit_lotsizingent = t; }
            else if (par == "time_limit_lotsizingdyn") { p.time_limit_lotsizingdyn = t; }
            else if (par == "time_limit_compact") { p.time_limit_compact = t; }
            else if (par == "time_limit_relax") { p.time_limit_relax = t; }
        }
    }
    else{
        cout << "erreur lecture fichier" << endl;
    }
}
