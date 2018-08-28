#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(){
    int nbOptGenCol = 0;
    int nbOptCompact = 0;
    float gapMoyGenCol = 0.0;
    float gapMoyCompact = 0.0;
    float tpsMoyGenCol = 0.0;
    float tpsMoyCompact = 0.0;
    for(int i=1; i<=288; ++i){
        string file = "../Log_gencolnostock/slurmout_"+to_string(i)+".out";
        ifstream fichier(file, ios::in); //ouverture du fichier
        if(fichier){  // si l'ouverture fonctionne
            string statutGenCol, statutCompact, ligneGenCol, ligneCompact;
            //getline(fichier,ligneGenCol);
            float gapGenCol, gapCompact, tpsGenCol, tpsCompact, sol;
            fichier >> sol >> tpsGenCol >> statutGenCol >> gapGenCol;
            getline(fichier,ligneCompact);
            fichier >> sol >> tpsCompact >> statutCompact >> gapCompact;
            if(statutGenCol == "optimal") nbOptGenCol++;
            if(statutCompact == "optimal") nbOptCompact++;
            else gapMoyCompact += gapCompact;
            gapMoyGenCol += gapGenCol;
            tpsMoyCompact += tpsCompact;
            tpsMoyGenCol += tpsGenCol;
        }
        else cout << "erreur ouverture fichier " << file << endl;
    }
    cout << "gapcompact " << gapMoyCompact << endl;
    cout << nbOptGenCol << " " << tpsMoyGenCol/288 << " " << gapMoyGenCol*100/288 << " // GENCOL : nbOpt, tpsMoy, gapMoy" << endl;
    cout << nbOptCompact << " " << tpsMoyCompact/288 << " " << gapMoyCompact*100/288 << " // COMPACT : nbOpt, tpsMoy, gapMoy" << endl;
    return 0;
}