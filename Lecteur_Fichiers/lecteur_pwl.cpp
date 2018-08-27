#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include "../struct.h"
#include "lecteur_pwl.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;

int inf = numeric_limits<int>::max();

void lecteur_pwl(string file, data &d){
	//d.cardT = cardT;
	//char temp;
	//vector< vector<int> > conflitstemp;
	


	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		//string ligne;
		//int t;
		float bp;
		float valbp;
		float pent;
		int nb_bp;
		fichier >> nb_bp;
		vector<int> nbp (d.cardT, nb_bp);
		/*vector<vector<float> > pente (d.cardT, vector<float> ());
		vector<vector<float> > valbpt (d.cardT, vector<float> ());
		vector<vector<float> > bpt (d.cardT, vector<float> ());*/
		vector<float> pente;
		vector<float> valbpt;
		vector<float> bpt;

		for(int i=0; i<nb_bp; ++i){ // tant que l'on peut lire une ligne
			fichier >> bp >> valbp >> pent;
			//cout << bp << " " << valbp << " " << pent << endl;
			pente.push_back(pent);
			valbpt.push_back(valbp);
			bpt.push_back(bp);
		}
		vector<vector<float> > vpente (d.cardT, pente);
		vector<vector<float> > vvalbpt (d.cardT, valbpt);
		vector<vector<float> > vbpt (d.cardT, bpt);

		//cout << vbpt[0][0] << " " << vbpt[1][0] << endl;


		d.pente = vpente;
		d.valbpt = vvalbpt;
		d.bpt = vbpt;
		d.nb_bp = nbp;

		
		for(int i=0; i<d.cardT; ++i){
			d.bpt[i].push_back(inf);
		}
		
		fichier.close();

		cout << "lecture pwl ok" << endl;
	}
	else{
		cout << "erreur lecture fichier pwl" << endl;
	}
}
