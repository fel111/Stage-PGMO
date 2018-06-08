#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
#include <limits>
//#include <math.h>
//#include <algorithm>
//#include <chrono>
//#include <deque>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "../struct.h"
#include "struct_gencol.h"
#include "../Lecteur_Fichiers/lecteur_taches.h"
#include "../Lecteur_Fichiers/lecteur_param.h"
#include "probScip.cpp"
#include "pricer.h"

using namespace std;
float infini = numeric_limits<float>::max();



void firstSol(structGenCol & sGC){
    int cptId = 0;
    
    for(int t=sGC.d.releaseDateMin; t<sGC.d.cardT; ++t){
        vector<int> taskList;
        feasibleSet l;
        for(int i=0; i<sGC.d.cardJ; ++i){
            if((sGC.d.rj[i]<=t) && (t<sGC.d.rj[i]+sGC.d.pj[i])){
                taskList.push_back(i);
                if(taskList.size() == 1){
                    l.energyDemand = sGC.d.Djk[i][0];
                    l.deadLine = sGC.d.dj[i];
                    l.releaseTime = sGC.d.rj[i];
                }
                else{
                    l.energyDemand += sGC.d.Djk[i][0];
                    if(l.deadLine > sGC.d.dj[i]) l.deadLine = sGC.d.dj[i];
                    if(l.releaseTime < sGC.d.rj[i]) l.releaseTime = sGC.d.rj[i];
                }
            }
        }
        if(taskList.size() > 0){
            l.id = cptId;
            l.timeGen = t; 
            l.tasksList = taskList;
            if(checkSet(l,sGC)){
                sGC.L.push_back(l);
                sGC.cardL = sGC.L.size();
                addSetK_l(l,sGC);
                addA_il(l,sGC);
                addL_t(l,sGC);
                cptId++;
            }else cout << "l deja present " << l.id << endl;
        }
	}

}


int initData(structGenCol & sGC){
    data d;
    param p;

    d.cardR = 2;
	d.s0 = 0;
	d.cardM = 1;
    string param = "param1.txt";
    string instance = "inst_2";
	if(lecteur_param("Param/"+param,p) == 0) return 0;
	
	p.qmax = 20;
	p.qmin = 0;
	p.qinit = 0;
	
	lecteur_taches_EnergSchedInst("Donnees/dataSchedulingInstances_newname/"+instance,d);

	initCap(d,p);

    for(int i=0;i<d.cardT;++i){
        //vector<float> temp_pente = {1.0,0.0,0.5,0.0,2.0,0.0,1.5};
		/*vector<float> temp_bpt = {0.0,4.0,4.0,8.0,8.0,100.0,100.0,infini};
		vector<float> temp_valbpt = {0.0,4.0,4.0,6.0,6.0,190.0,190.0};
		vector<float> temp_pente = {1.0,0.5,2.0,1.5};*/
		//vector<float> temp_bpt = {0.0,4.0,8.0,100.0,infini};
		//vector<float> temp_valbpt = {0.0,4.0,6.0,190.0};
		vector<float> temp_bpt = {0.0,5.0,5.0,1000.0};
		vector<float> temp_valbpt = {0.0,10.0,10.0,3000.0};
		vector<float> temp_pente = {2.0,0.5,2.0};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(4);

        //initialisation L_t
        vector<int> lt_temp;
        sGC.L_t.push_back(lt_temp);
	}
	vector<vector<float> > ress;
	for(int i=0;i<d.cardM;++i){
		//vector<float> r = {(d.cardJ*4.0/d.cardM),(d.cardJ*4000.0/d.cardM)}; // s'assure qu'il y ai assez de ram/cpu
		vector<float> r = {d.cardJ, d.cardJ};
		//cout << r[0] << " " << r[1] << endl;
		ress.push_back(r);
	}
	d.Ckr = ress;
	vector<float> Dk (d.cardM, 0.0); // ne converge pas avec != 0
	d.Dk = Dk;

    sGC.d = d;
    sGC.p = p;

    //initialisation a_il
    for(int j=0; j<d.cardJ; ++j){
        vector<int> ailtemp;
        sGC.a_il.push_back(ailtemp);
    }

}


SCIP_RETCODE ColGen_Scip(structGenCol & sGC){
	int tmp;
    SCIP_RETCODE status;
    SCIP_STATUS solstatus;
    initData(sGC);
    //cout <<"initData OK ! "<<endl;
	
	firstSol(sGC);
	//cout <<"firstSol OK ! "<<endl;
	
	//affK_l(sGC);
	//affL_t(sGC);
	affAllSet(sGC);
	
    // Load and build model
    SCIP_CALL (status = Load_Original_Model(sGC));
    
    // Set dedicated parameters
    //SCIP_CALL (status = SetColGenParameters_Scip((*pU).scip, (*pU).Params));
    
    
    // activate the pricer
    SCIP_CALL( SCIPactivatePricer(sGC.scip, SCIPfindPricer(sGC.scip, "generatetestnewobjects") ));
    
    // SOLVE
    //cout << "Start solve..."<<endl;
    SCIP_CALL( status = SCIPsolve(sGC.scip) );


    sGC.sol = SCIPgetBestSol(sGC.scip);

    cout << "verif : "<<verifSol(sGC) << endl; 

    for(int j=0; j<sGC.d.cardJ; ++j){
        for(int t=sGC.d.rj[j]; t<sGC.d.dj[j]; ++t){
            if(SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][t]) == 1) cout <<"tache "<<j<<" réalisée au temps "<<t<<endl;
        }
    }
    //cout<<"... End solve"<<endl;
    /*
    printf("\nstart solution export ...\n");
    // Get the current scip solving time in seconds.
    (*pU).ScipResolutionTime = SCIPgetSolvingTime ((*pU).scip);
    
    // Get solution
    (*pU).Solution = SCIPgetBestSol((*pU).scip);
    assert( (*pU).Solution != NULL);
    
    // Print results in log and summary files
    solstatus = SCIPgetStatus((*pU).scip) ;
    SCIP_CALL(status = PrintColGenResults (pU, solstatus));
    printf("\n...end solution export\n");
    
    // Free memory
    assert(tmp == 0);
    */
    return status;
}




int main(){
    structGenCol sGC;

    ColGen_Scip(sGC);

    return 0;
}