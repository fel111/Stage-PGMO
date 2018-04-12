#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void fichier_sol(int cardT, int nb_bp, int bornesup, int stockmax, int nbfich, float tpsDyn, float solDyn, float tpsPLNE, float solPLNE, float tpsMIP, float solMIP, float binfPLNE, float binfMIP, string statusPLNE, string statusMIP){
    string fichier = "../Log_sol/"+to_string(cardT)+"_"+to_string(nb_bp)+"_"+to_string(bornesup)+"_"+to_string(stockmax)+"_"+to_string(nbfich);

    ofstream f(fichier, ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier

    if(f){  // si l'ouverture a réussi
        f << solDyn << " " << tpsDyn << " //temps et sol LotSizDyn" << endl;
        f << solPLNE << " " << tpsPLNE <<" "<<binfPLNE<<" "<<statusPLNE <<" //temps, sol, binf et status PLNE CPLEX" << endl;
        f << solMIP << " " << tpsMIP << " "<<binfMIP<<" "<<statusMIP<< " //temps, sol, binf et status MIP CPLEX" << endl;
        //f << (solDyn - solMIP)*100/solMIP << " //ecart relatif entre solDyn et solOpt (%)"<<endl;
        		
        f.close();  // on referme le fichier
    }

    else cerr << "Erreur à l'ouverture !" << endl;


}
