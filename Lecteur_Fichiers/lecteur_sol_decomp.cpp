#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
//#include "struct.h"
//#include <algorithm>
//#include <chrono>
//#include <random>
//#include <list>
using namespace std;


void lecteur_sol(int param){

    int nbfichiersok=0;
	float tpscpxmoy=0;
    float tpscpxmin;
    float tpscpxmax;
    float tpsordomoy=0;
    float tpsordomin;
    float tpsordomax;
    float tpslsmoy=0;
    float tpslsmin;
    float tpslsmax;
    float bornesupmoycpx=0;
    float bornesupmoydec=0;
    float nboptcpx=0;
    float nboptdec=0;
    bool first = true;

    for(int i=1; i<=288; ++i){


        ifstream fichier("../slurmouts_sansrelax/slurmout_"+to_string(i)+"_"+to_string(param)+".out", ios::in); //ouverture du fichier
        if(fichier){  // si l'ouverture fonctionne
            int nbligne=0;
            string ligne;
            


            while(getline(fichier,ligne)) ++nbligne;

            if(nbligne==8){
                ++nbfichiersok;
                string inut;
                string tpscpx;
                string statuscpx;
                string bsupcpx;
                string tpsdec;
                string tpsordo;
                string tpsls;
                string bsupdec;
                //reinit pour getline
                fichier.clear();
                fichier.seekg(0, ios::beg);

                getline(fichier,ligne);
                getline(fichier,ligne);
                istringstream iss (ligne);
                iss >> inut >> tpscpx >> statuscpx >> bsupcpx;
                getline(fichier,ligne);
                getline(fichier,ligne);
                getline(fichier,ligne);
                getline(fichier,ligne);
                istringstream iss2 (ligne);
                iss2 >> inut >> inut >> tpsdec >> tpsordo >> tpsls >> bsupdec;
        
                tpscpx.erase(0,4);
                tpsdec.erase(0,10);
                tpsordo.erase(0,8);
                tpsls.erase(0,6);
                bsupcpx.erase(0,9);
                bsupdec.erase(0,9);

                if(bsupdec != "-1"){

                    //cout << tpscpx << " " << tpsdec << " " << tpsordo << " " << tpsls << " " << bsupcpx << " " << bsupdec << endl;

                    if(statuscpx == "statut=Optimal"){ 
                        nboptcpx++; 
                        if(bsupdec==bsupcpx) nboptdec++;
                    }
                    if(first){
                        tpscpxmax=stof(tpscpx,0);
                        tpscpxmin=stof(tpscpx,0);
                        tpsordomax=stof(tpsordo,0);
                        tpsordomin=stof(tpsordo,0);
                        tpslsmax=stof(tpsls,0);
                        tpslsmin=stof(tpsls,0);
                    }
                    else{
                        if(stof(tpscpx,0)>tpscpxmax) tpscpxmax=stof(tpscpx,0);
                        if(stof(tpscpx,0)<tpscpxmin) tpscpxmin=stof(tpscpx,0);
                        if(stof(tpsordo,0)>tpsordomax) tpsordomax=stof(tpsordo,0);
                        if(stof(tpsordo,0)<tpsordomin) tpsordomax=stof(tpsordo,0);
                        if(stof(tpsls,0)>tpslsmax) tpslsmax=stof(tpsls,0);
                        if(stof(tpsls,0)<tpslsmin) tpslsmin=stof(tpsls,0);
                    }
                    tpscpxmoy+=stof(tpscpx,0);
                    tpslsmoy+=stof(tpsls,0);
                    tpsordomoy+=stof(tpsordo,0);
                    bornesupmoycpx+=stof(bsupcpx,0);
                    bornesupmoydec+=stof(bsupdec,0);
                    

                    first=false;
                }
                else nbfichiersok--;
            }
            fichier.close();
        }
        else{
            cout << "erreur lecture fichier " << i << endl;
        }
    }

    tpscpxmoy = tpscpxmoy/nbfichiersok;
    tpslsmoy = tpslsmoy/nbfichiersok;
    tpsordomoy = tpsordomoy/nbfichiersok;
    float tpsdecmoy = tpslsmoy+tpsordomoy;

    cout << "tpsCPLEXmoy = " << tpscpxmoy << "   tpsDECmoy = " << tpsdecmoy << endl;
    cout << "nb opt CPX = " << nboptcpx << "   nb opt dec = " << nboptdec << endl;
    cout << "nb fichiers ok = " << nbfichiersok << endl;
    cout << "ratio bsupCPX/bsupDEC = " << (bornesupmoydec-bornesupmoycpx)/bornesupmoycpx*100 << endl;
}

int main(int argc, char *argv[]){
    
    lecteur_sol(10);
}