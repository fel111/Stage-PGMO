#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>

using namespace std;

int main (int argc, char *argv[]){
	if(argc != 3){ 
		cout << " Paramètres demandés : nombre périodes, borne demande" <<endl;
		return 0;
	}
	int cardT = atoi(argv[1]);
	int max = atoi(argv[2]);

	
	string fichier = "demande_"+to_string(cardT)+"_"+to_string(max);

        ofstream f(fichier, ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier

        if(f){  // si l'ouverture a réussi

		vector<int> dt;
		random_device rd;
    		mt19937 mt(rd());
    		uniform_int_distribution<int> dist(1, max);
		for(int i=0; i<cardT; ++i){
			dt.push_back(dist(mt));
			f << i << " " << dt[i] << endl;
		}	



        	f.close();  // on referme le fichier
	}

        else cerr << "Erreur à l'ouverture !" << endl;

 

        return 0;

}
