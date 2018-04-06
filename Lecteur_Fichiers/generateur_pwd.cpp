#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>

using namespace std;

void generateur_pwd(int cardT, int nb_bp, int pentemax, int bornesup){
	
	string fichier = "pwd_"+to_string(cardT)+"_"+to_string(nb_bp)+"_"+to_string(bornesup);

        ofstream f(fichier, ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier

        if(f){  // si l'ouverture a réussi
		f << cardT << " " << nb_bp << endl;
		vector<float> pente;
		vector<int> bpt;
		vector<float> valbpt;
		random_device rd;
    		mt19937 mt(rd());
    		uniform_real_distribution<float> dist(0.5, (float) pentemax);
		for(int i=0; i<nb_bp; ++i){
			pente.push_back(floor(dist(mt) * 10.) / 10.);
		}
		bpt.push_back(0);
		valbpt.push_back(0.0);
		int dx = bornesup / nb_bp;
		for(int i=1; i<nb_bp; ++i){
			bpt.push_back(dx*i);
			valbpt.push_back((bpt[i]-bpt[i-1])*pente[i-1]+valbpt[i-1]);
		}
		for(int k=0; k<cardT; ++k){
			for(int i=0; i<nb_bp; ++i){
				f << k << " " << bpt[i] << " " << valbpt[i] << " " << pente[i] << endl;
			}	
		}		
		



        	f.close();  // on referme le fichier
	}

        else cerr << "Erreur à l'ouverture !" << endl;


}
