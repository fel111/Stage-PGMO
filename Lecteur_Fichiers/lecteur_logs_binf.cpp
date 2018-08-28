#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(){
    int nbOpt = 0;
    int nbBestG = 0;
    int nbEq = 0;
    float gap = 0.0;
    for(int i=1; i<=288; ++i){
        if((i!=126)&&(i!=132)){
            string file = "../Log_binf_gencol_compact_600/slurmout_"+to_string(i)+".out";
            ifstream fichier(file, ios::in); //ouverture du fichier
            if(fichier){  // si l'ouverture fonctionne
                string binfgencol,ligneCompact;
                float binfG, binfC;
                fichier >> binfgencol;
                //cout << binfgencol << " " << file<< endl;

                if(binfgencol != "-1e+20"){
                    getline(fichier,ligneCompact);
                    fichier >> binfC;
                    
                    binfG = stof(binfgencol);
                    if (binfG > binfC) nbBestG++;
                    if (binfG == binfC) nbEq++;
                // cout << binfC << " " << binfG << endl;
                    ++nbOpt;
                    gap += (binfG-binfC)/binfC;
                }
            }
            else cout << "erreur ouverture fichier " << file << endl;
        }
    }
    cout << "gapmoyen " << (gap*100)/nbOpt << endl;
    cout << "nbOPt : " << nbOpt << endl;
    cout << "nbBestG : " << nbBestG << "   nbEq : " << nbEq << endl;
    return 0;
}