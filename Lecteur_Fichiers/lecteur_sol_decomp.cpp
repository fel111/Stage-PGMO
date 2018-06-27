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


void lecteur_sol(int taille){

    
    //float tpscpxmin;
    //float tpscpxmax;
    //float tpsordomoy=0;
    //float tpsordomin;
    //float tpsordomax;
   // float tpslsmoy=0;
   // float tpslsmin;
    //float tpslsmax;
    //float bornesupmoycpx=0;
    //float bornesupmoydec=0;
    //float nboptcpx=0;
    //float nboptdec=0;
    bool first = true;
    vector<int> capacite = {5,10,15,25};


    for(const auto& cap : capacite){

        int nbfichiersok3=0;
        int nbfichiersok5=0;
        int nbfichiersok10=0;
        int nbconv3=0;
        int nbconv5=0;
        int nbconv10=0;
        int nbboucleconv3=0;
        int nbboucleconv5=0;
        int nbboucleconv10=0;
        int nboptcpx=0;
        float gapmoy3=0.0;
        float gapmoy5=0.0;
        float gapmoy10=0.0;
        float tpsmoy3=0.0;
        float tpsmoy5=0.0;
        float tpsmoy10=0.0;
        float tpsmoyordo3=0.0;
        float tpsmoyls3=0.0;
        float tpsmoyordo5=0.0;
        float tpsmoyls5=0.0;
        float tpsmoyordo10=0.0;
        float tpsmoyls10=0.0;
        float tpsmoycpx=0.0;
        bool firstcpx=true;

        vector<int> fichierparam;
        if(cap==5) fichierparam = {1,5,9};
        if(cap==10) fichierparam = {2,6,10};
        if(cap==15) fichierparam = {3,7,11};
        if(cap==25) fichierparam = {4,8,12};

        for(const auto& f : fichierparam){
            
            for(int i=1; i<=288; ++i){
                cout << "i, f : " <<i<< " "<<f<<endl;
                ifstream fichier("../slurmouts_sansrelax/slurmout_"+to_string(i)+"_"+to_string(f)+".out", ios::in); //ouverture du fichier
                if(fichier){  // si l'ouverture fonctionne
                    int nbligne=0;
                    string ligne;
                


                    while(getline(fichier,ligne)) ++nbligne;

                    if(nbligne==8){
                        // ·++nbfichiersok;
                        string inut;
                        string tpscpx;
                        string statuscpx;
                        string bsupcpx;
                        string tpsdec;
                        string tpsordo;
                        string tpsls;
                        string bsupdec;
                        string statutdec;
                        string nbboucledec;
                        string tpsdecpl;
                        string tpsordopl;
                        string tpslspl;
                        string bsupdecpl;
                        string statutdecpl;
                        string nbboucledecpl;
                        string nbtaches;
                        //reinit pour getline
                        fichier.clear();
                        fichier.seekg(0, ios::beg);

                        getline(fichier,ligne);
                        istringstream iss4 (ligne);
                        iss4 >> inut >> inut >> inut >> nbtaches;
                        getline(fichier,ligne);
                        istringstream iss (ligne);
                        iss >> inut >> tpscpx >> statuscpx >> bsupcpx;
                        getline(fichier,ligne);
                        getline(fichier,ligne);
                        istringstream iss2 (ligne);
                        iss2 >> inut >> inut >> inut >> tpsdecpl >> tpsordopl >> tpslspl >> bsupdecpl >> nbboucledecpl >> statutdecpl;
                        getline(fichier,ligne);
                        getline(fichier,ligne);
                        istringstream iss3 (ligne);
                        iss3 >> inut >> inut >> tpsdec >> tpsordo >> tpsls >> bsupdec >> nbboucledec >> statutdec;
                
                        tpscpx.erase(0,4);
                        tpsdec.erase(0,10);
                        tpsordo.erase(0,8);
                        tpsls.erase(0,6);
                        bsupcpx.erase(0,9);
                        bsupdec.erase(0,9);
                        tpsdecpl.erase(0,10);
                        tpsordopl.erase(0,8);
                        tpslspl.erase(0,6);
                        bsupdecpl.erase(0,9);
                        nbboucledecpl.erase(0,9);
                        nbboucledec.erase(0,9);

                        if(nbtaches==("nbTaches="+to_string(taille))){
                            if(firstcpx){
                                if(statuscpx == "statut=Optimal") nboptcpx++;
                                tpsmoycpx+=stof(tpscpx,0);
                            }

                            if(bsupdec != "-1"){
                                //++nbfichiersok;
                                if(f%4 == 0){
                                    if(statutdec=="arret=CONV"){
                                        nbconv3++;
                                        nbboucleconv3 += stoi(nbboucledec,0);
                                    }
                                    ++nbfichiersok3;
                                    tpsmoy3+=stof(tpsdec,0);
                                    tpsmoyordo3+=stof(tpsordo,0);
                                    tpsmoyls3+=stof(tpsls,0);
                                    gapmoy3+=((stof(bsupdec,0)-stof(bsupcpx,0))/stof(bsupcpx,0));
                                }
                                else if(f%4 == 1){
                                    if(statutdec=="arret=CONV"){
                                        nbconv5++;
                                        nbboucleconv5+=stoi(nbboucledec,0);
                                    }
                                    ++nbfichiersok5;
                                    tpsmoy5+=stof(tpsdec,0);
                                    tpsmoyordo5+=stof(tpsordo,0);
                                    tpsmoyls5+=stof(tpsls,0);
                                    gapmoy5+=((stof(bsupdec,0)-stof(bsupcpx,0))/stof(bsupcpx,0));
                                }
                                else if(f%4 == 2){
                                    if(statutdec=="arret=CONV"){
                                        nbconv10++;
                                        nbboucleconv10+=stoi(nbboucledec,0);
                                    }
                                    ++nbfichiersok10;
                                    tpsmoy10+=stof(tpsdec,0);
                                    tpsmoyordo10+=stof(tpsordo,0);
                                    tpsmoyls10+=stof(tpsls,0);
                                    gapmoy10+=((stof(bsupdec,0)-stof(bsupcpx,0))/stof(bsupcpx,0));
                                }
                            }
                        }
                    }
                    fichier.close();
                }
                else{
                    cout << "erreur lecture fichier " << i << endl;
                }
            }
            firstcpx=false;
        }

        tpsmoycpx = tpsmoycpx/144;
        tpsmoy3 = tpsmoy3/nbfichiersok3;
        tpsmoy5 = tpsmoy5/nbfichiersok5;
        tpsmoy10 = tpsmoy10/nbfichiersok10;
        tpsmoyls3 = tpsmoyls3/nbfichiersok3;
        tpsmoyls5 = tpsmoyls5/nbfichiersok5;
        tpsmoyls10 = tpsmoyls10/nbfichiersok10;
        tpsmoyordo3 = tpsmoyordo3/nbfichiersok3;
        tpsmoyordo5 = tpsmoyordo5/nbfichiersok5;
        tpsmoyordo10 = tpsmoyordo10/nbfichiersok10;
        gapmoy3 = gapmoy3/nbfichiersok3;
        gapmoy5 = gapmoy5/nbfichiersok5;
        gapmoy10 = gapmoy10/nbfichiersok10;

        cout << cap << " : | " << nbconv3 << " " <<  tpsmoy3 << " " << gapmoy3 << " | " << nbconv5 << " " <<  tpsmoy5 << " " << gapmoy5 << " | " << nbconv10 << " " <<  tpsmoy10 << " " << gapmoy10 << " | " << tpsmoycpx << " " << nboptcpx << endl;
    }
}

int main(int argc, char *argv[]){
    
    lecteur_sol(30);
}