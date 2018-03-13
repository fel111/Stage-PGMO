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
		int tmp = -1;
		int t;
		int bp;
		int valbp;
		float pentepr;
		float penteder;
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t;
			if(t != tmp){				
				iss >> pentepr >> penteder;
				pente[0].push_back(pentepr);
				pente[1].push_back(penteder);
				tmp = t;	
			}
			else{
				iss >> bp >> valbp;
				bpt[t].push_back(bp);
				valbpt[t].push_back(valbp);
			}
		}
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}
