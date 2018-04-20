#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include "data_struct.h"
#include "lecteur_pwd.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;

int inf = numeric_limits<int>::max();

void lecteur_pwd(string file, data &d){
	//d.cardT = cardT;
	//char temp;
	//vector< vector<int> > conflitstemp;
	


	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		float bp;
		float valbp;
		float pent;
		int nb_bp;
		getline(fichier,ligne);
		istringstream iss (ligne);
		iss >> d.cardT >> nb_bp;
		vector<int> nbp (d.cardT, nb_bp);
		vector<vector<float> > pente (d.cardT, vector<float> ());
		vector<vector<float> > valbpt (d.cardT, vector<float> ());
		vector<vector<float> > bpt (d.cardT, vector<float> ());


		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> bp >> valbp >> pent;
			bpt[t].push_back(bp);
			valbpt[t].push_back(valbp);
			pente[t].push_back(pent);
		}
		d.pente = pente;
		d.valbpt = valbpt;
		d.bpt = bpt;
		d.nb_bp = nbp;

		fichier.close();

		for(int i=0; i<d.cardT; ++i){
			d.bpt[i].push_back(inf);
		}

	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}
