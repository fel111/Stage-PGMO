#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <limits>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <deque>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "data_struct.h"
#include "Ordo_compact/ordo_comp.h"
#include "Lot Sizing/Lot Sizing Dynamique/lotsizdyn2.h"
//#include <ilcplex/cplex.h>
//#define SCIP_DEBUG
using namespace std;
int infini_int = numeric_limits<int>::max();


int main(){
	data d;
	// init test

	d.cardJ = 4;
	d.cardM = 2;
	d.cardR = 2;
	d.s0 = 10;
	d.cardT = 10; 
	d.Q = 20;
	for(int i=0;i<d.cardT;++i){
		vector<float> temp_pente = {3.0,2.0,1.0};
		vector<int> temp_bpt = {0,10,20,infini_int};
		vector<float> temp_valbpt = {0.0,30.0,50.0};
		d.pente.push_back(temp_pente);
		d.bpt.push_back(temp_bpt);
		d.valbpt.push_back(temp_valbpt);
		d.nb_bp.push_back(3);
	}

	d.pj = {4,2,8,5};
	d.cjr = {{20.0,1.0},{15.0,1.0},{35.0,1.0},{10.0,1.0}};
	d.Ckr = {{300.0,4.0},{500.0,2.0}};
	d.Dk = {10.0,15.0};
	d.Djk = {{10.0,20.3},{12.1,30.2},{12.4,5.2},{17.1,54.6}};

	d.dt = ordo(d);
	cout << "ORDO REALISE" << endl;
	for(int i=0; i<d.cardT; ++i) cout << d.dt[i] << "," << endl;
	


	return 0;
}








