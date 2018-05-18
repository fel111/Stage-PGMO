#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include "struct.h"
#include "lecteur_taches.h"
//#include <algorithm>
//#include <chrono>
#include <random>
//#include <list>
using namespace std;


void lecteur_taches(string file, data &d, bool fenetre){ //bool precise si l'on souhaite fenetre de temps ou non

	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		float conso,cpu,ram;
        int duree, releasedate, duedate;
        int cardJ = d.cardJ;
        vector<vector<float> > cjr;
        vector<vector<float> > Djk;
        vector<int> pj, rj, dj;
		while((getline(fichier,ligne))&&(cardJ>0)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> duree >> releasedate >> duedate >> conso >> cpu >> ram;
            vector<float> ress;
            ress.push_back(cpu);
            ress.push_back(ram);
            cjr.push_back(ress);
			vector<float> cons (d.cardM,conso);
            Djk.push_back(cons);
            pj.push_back(duree);
            if(fenetre){
                rj.push_back(releasedate);
                dj.push_back(duedate);
            }
            else{
                rj.push_back(0);
                dj.push_back(d.cardT-1);
            }
            --cardJ;
		}
		d.Djk = Djk;
        d.cjr = cjr;
        d.pj = pj;
        d.rj = rj;
        d.dj = dj;
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}


void lecteur_taches_EnergSchedInst(string file, data &d){ //bool precise si l'on souhaite fenetre de temps ou non

	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
        float conso;
        float consoTot = 0.0;
        int duree, releasedate, duedate, cardJ;
        int cardT = 0;     
        vector<vector<float> > Djk;
        vector<int> pj, rj, dj;
        getline(fichier,ligne);
        istringstream iss1 (ligne);
        iss1 >> cardJ;
        vector<vector<float> > cjr (cardJ, vector<float> (2,1.0));
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss2 (ligne);
			iss2 >> conso >> duree >> releasedate >> duedate;
            if(duedate > cardT) cardT = duedate;
			vector<float> cons (d.cardM,conso);
            Djk.push_back(cons);
            consoTot += conso;
            pj.push_back(duree);
            rj.push_back(releasedate);
            dj.push_back(duedate);
		}
		d.Djk = Djk;
        d.cjr = cjr;
        d.pj = pj;
        d.rj = rj;
        d.dj = dj;
        d.cardJ = cardJ;
        d.cardT = cardT+1;
        d.consoTot = consoTot;
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}



void generateur_taches(string file, int cardJ, int cardT, float consomax){

    ofstream f(file, ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier

    if(f){  // si l'ouverture a réussi
        random_device rd;
        mt19937 mt(rd());
        uniform_real_distribution<float> distfloat(1, consomax);
        uniform_int_distribution<int> distduree(1, cardT-1);
        uniform_int_distribution<int> distcpu(1, 4);
        uniform_int_distribution<int> distram(1, 4000);
        for(int i=0; i<cardJ; ++i){
            int duree = distduree(mt);
            cout << " duree "<<duree << endl;
            uniform_int_distribution<int> distrj(1, cardT-duree);
            int rj = distrj(mt);
            cout << rj << endl;
            uniform_int_distribution<int> distdj(rj+duree, cardT);
            int dj = distdj(mt);
            cout<< dj<< endl;
            f << duree << " " << rj << " " << dj << " " << distfloat(mt) << " " << distcpu(mt) << " "<< distram(mt) << endl;	
        }


        f.close();  // on referme le fichier
    }
    else cerr << "Erreur à l'ouverture !" << endl;

}

/*int main(){
    generateur_taches("../Donnees/taches_120",200,120,2000);

    return 0;
}*/
