#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <random>
#include <cmath>
#include "data_struct.h"

using namespace std;

/*float roundd(double var, int nbdec)
{
    // 37.66666 * 100 =3766.66
    // 3766.66 + .5 =37.6716    for rounding off value
    // then type cast to int so value is 3766
    // then divided by 100 so the value converted into 37.66
    float value = (int)(var * pow(10,nbdec) + .5);
    return (float)value / pow(10,nbdec);
}*/

void generateur_demande(int cardT, int max){
	/*if(argc != 3){ 
		cout << " Paramètres demandés : nombre périodes, borne demande" <<endl;
		return 0;
	}
	int cardT = atoi(argv[1]);
	int max = atoi(argv[2]);*/
	for(int j=1; j<11; ++j){
	
		string fichier = "../Donnees/demande_"+to_string(cardT)+"_"+to_string(max)+"_"+to_string(j);

			ofstream f(fichier, ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier

			if(f){  // si l'ouverture a réussi
			f << cardT << endl;
			vector<float> dt;
			random_device rd;
			mt19937 mt(rd());
			uniform_real_distribution<float> dist(1, max);
			for(int i=0; i<cardT; ++i){
				dt.push_back(dist(mt));
				f << i << " " << roundd(dt[i],1) << endl;
			}	



				f.close();  // on referme le fichier
		}

			else cerr << "Erreur à l'ouverture !" << endl;
	}

}

