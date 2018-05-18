#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "struct.h"
#include "lecteur_demande.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;


void lecteur_demande(string file, data &d){
	//d.cardT = cardT;
	//char temp;
	//vector< vector<int> > conflitstemp;
	



	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		float dem;
        getline(fichier,ligne);
        istringstream iss (ligne);
        iss >> t;
		d.cardT = t;
		vector<float> dt (d.cardT, 0.0);
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> dem;
			dt[t] = dem;
		}
		d.dt = dt;
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}
