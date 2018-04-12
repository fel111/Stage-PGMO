#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
//#include "data_struct.h"
#include "lecteur_sol.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;


void lecteur_sol(int cardT, int nbbp,int bornesup, int stockmax,string fichsol){
	//d.cardT = cardT;
	//char temp;
	//vector< vector<int> > conflitstemp;
	float tpsmoyDyn = 0.0;
    float tpsmoyEnt = 0.0;
    float tpsmoyMip = 0.0;
    float ecartmoy = 0.0;
    float tpsmaxDyn = 0.0;
    float tpsmaxEnt = 0.0;
    float tpsmaxMip = 0.0; 
    float tpsminDyn = 1000000.0;
    float tpsminEnt = 1000000.0;
    float tpsminMip = 1000000.0;
    float bsupEnt = 0.0;
    float binfEnt = 0.0;
    float bsupMip = 0.0;
    float binfMip = 0.0;
    float bsupDyn = 0.0;
    int nbDynFail = 0;
    int nbEntFail = 0;
    int nbMipFail = 0;
    int nbFailOpen = 0;
    string file = "../Log_sol/"+to_string(cardT)+"_"+to_string(nbbp)+"_"+to_string(bornesup)+"_"+to_string(stockmax);
    for(int i=1; i<11; ++i){


        ifstream fichier(file+"_"+to_string(i), ios::in); //ouverture du fichier
        if(fichier){  // si l'ouverture fonctionne
            string ligne;
            float soldyn,tpsdyn,solent,tpsent,solmip,tpsmip,binfm,binfe;


            getline(fichier,ligne);
            istringstream iss (ligne);
            iss >> soldyn >> tpsdyn;
            
            //cout << "soldyn : " << soldyn << endl;
            if(tpsdyn < 300.0){
                bsupDyn += soldyn;
                tpsmoyDyn += tpsdyn;
                if(tpsdyn > tpsmaxDyn) tpsmaxDyn = tpsdyn;
                if(tpsdyn < tpsminDyn) tpsminDyn = tpsdyn;
            }
            else ++nbDynFail;

            getline(fichier,ligne);
            istringstream isss (ligne);
            isss >> solent >> tpsent >> binfe;
            
            //cout << "soldyn : " << soldyn << endl;
            tpsmoyEnt += tpsent;
            bsupEnt += solent;
            binfEnt += binfe;
            if(tpsent > tpsmaxEnt) tpsmaxEnt = tpsent;
            if(tpsent < tpsminEnt) tpsminEnt = tpsent;

            getline(fichier,ligne);
            istringstream issss (ligne);
            issss >> solmip >> tpsmip >> binfm;
            //cout << "soldyn : " << soldyn << endl;
            tpsmoyMip += tpsmip;
            bsupMip += solmip;
            binfMip += binfm;
            if(tpsmip > tpsmaxMip) tpsmaxMip = tpsmip;
            if(tpsmip < tpsminMip) tpsminMip = tpsmip;

            if(tpsdyn < 300.0) ecartmoy += (soldyn - solmip)*100/solmip;
            //else ++nbDynFail;
            if(tpsent>300.0) ++nbEntFail;
            if(tpsmip>300.0) ++nbMipFail;

            
            fichier.close();
        }
        else{
            ++nbFailOpen;
            cout << "erreur lecture fichier " << i << endl;
        }
    }
    //cout << binfMip/10<< endl;
    if(nbFailOpen == 0){
    ofstream f(fichsol, ios::out | ios::app);
    //cout << "\multicolumn{1}{|c|}{"<<cardT <<"} & \multicolumn{1}{c|}{"<<nbbp <<"} &"<< stockmax<<"&"<<tpsmoyDyn/(10-nbDynFail) <<"&"<<tpsminDyn<<"&"<<tpsmaxDyn<<"&"<<10-nbDynFail <<"/10 &"<<tpsmoyEnt/10 <<"&"<<tpsminEnt <<"&"<<tpsmaxEnt <<"&"<<bsupEnt/10 <<"&"<<binfEnt/10 <<"&"<<10-nbEntFail <<"&"<<tpsmoyMip/10 <<"&"<<tpsminMip<<"&"<<tpsmaxMip <<"&"<<bsupMip/10 <<"&"<<binfMip/10 <<"&"<<10-nbMipFail <<"\\ \hline"<<endl;
    f<<cardT <<","<<nbbp <<","<< stockmax<<","<<tpsmoyDyn/(10-nbDynFail) <<","<< bsupDyn/(10-nbDynFail)<<","<<tpsminDyn<<","<<tpsmaxDyn<<","<<10-nbDynFail <<"/10 ,"<<tpsmoyEnt/10 <<","<<tpsminEnt <<","<<tpsmaxEnt <<","<<bsupEnt/10 <<","<<binfEnt/10 <<","<<10-nbEntFail <<"/10 ,"<<tpsmoyMip/10 <<","<<tpsminMip<<","<<tpsmaxMip <<","<<bsupMip/10 <<","<<binfMip/10 <<","<<10-nbMipFail<<"/10,"<<ecartmoy/(10-nbDynFail)<<endl;
    f.close();
    /*cout << "temps moyen Prog Dyn   : " << tpsmoyDyn/(10-nbDynFail) << "  min : "<< tpsminDyn << "  max : "<< tpsmaxDyn <<"   ("<<10-nbDynFail<<"/10 instances résolues avant 300s)" << endl;
    cout << "temps moyen CPLEX PLNE : " << tpsmoyEnt/10 << "  min : "<< tpsminEnt << "  max : "<< tpsmaxEnt <<"   ("<<10-nbEntFail<<"/10 instances résolues avant 300s)" << endl;
    cout << "temps moyen CPLEX MIP  : " << tpsmoyMip/10 << "  min : "<< tpsminMip << "  max : "<< tpsmaxMip <<"   ("<<10-nbMipFail<<"/10 instances résolues avant 300s)" << endl;
    cout << "ecart moyen Prog dyn et MIP : " << ecartmoy/(10-nbDynFail) << endl;*/
    }
}

int main(int argc, char *argv[]){
    /*string cardT = argv[1];
	string nb_bp = argv[2];
    string bornesup = argv[3];
	string stockmax = argv[4];*/
    vector<int> cardT = {20,100,500,1000};
    vector<int> stockmax = {100,1000,10000};
    vector<int> nbbp = {4,10};
    int bornesup = 2000;

    for(auto const& t : cardT){
        for(auto const& s : stockmax){
            for(auto const& bp : nbbp){
                lecteur_sol(t,bp,bornesup,s,"resultat4.csv");
            }

        }

    }


    
    //lecteur_sol(stoi(cardT),stoi(nb_bp),stoi(bornesup),stoi(stockmax),"resultat3.csv");
}