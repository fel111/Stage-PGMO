#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;


void lecture_fichier(string file, vector<vector<int> >& bpt, vector<vector<int> >& valbpt, vector<vector<float> >& pente){

	//char temp;
	//vector< vector<int> > conflitstemp;
	
	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		int bp;
		float valbp;
		float pente;
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> bp >> valbp >> pente;
			bpt[t].push_back(bp);
			valbpt[t].push_back(valbp);
			pente[t].push_back(pente);
		}
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}
