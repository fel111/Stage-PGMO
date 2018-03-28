#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>

using namespace std;

int main (int argc, char *argv[]){
	if(argc != 5){ 
		cout << " Paramètres demandés : nombre périodes, nombre breakpoints, pente max, bornesup" <<endl;
		return 0;
	}
	int cardT = atoi(argv[1]);
	int nb_bp = atoi(argv[2]);
	int pentemax = atoi(argv[3]);
	int bornesup = atoi(argv[4]);
	
	string fichier = "pwd_"+to_string(cardT)+"_"+to_string(nb_bp)+"_"+to_string(pentemax)+"_"+to_string(bornesup);

        ofstream f(fichier, ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier

        if(f){  // si l'ouverture a réussi
		vector<float> pente;
		vector<int> bpt;
		vector<float> valbpt;
		random_device rd;
    		mt19937 mt(rd());
    		uniform_real_distribution<float> dist(0.0, (float) pentemax);
		for(int i=0; i<nb_bp; ++i){
			pente.push_back(floor(dist(mt) * 100.) / 100.);
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

 

        return 0;

}
