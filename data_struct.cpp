#include <vector>
#include <math.h>
#include "data_struct.h"
using namespace std;


void modifPWL(data &d, vector<int> rt){
	for(int t=0; t<d.cardT; ++t){
		for(int i=0; i<d.nb_bp[t]; ++i){
			d.bpt[t][i] += rt[t];
			//if(rt[t] > 0.1) d.bpt[t][i] += rt[t];
			//else if(rt[t]<0) d.bpt[t][i] -= rt[t]; 
		}
	}
}

void modifPWL(data &d, vector<float> rt){
	for(int t=0; t<d.cardT; ++t){
		for(int i=0; i<d.nb_bp[t]; ++i){
			//if (rt[t] > 0.5) d.bpt[t][i] += rt[t];
			d.bpt[t][i] += rt[t];
			//if(rt[t] > 0) d.bpt[t][i] += rt[t];
			//else if(rt[t]<0) d.bpt[t][i] -= rt[t]; 
		}
	}
}

vector<int> dtToInt(data d, int choix){
	vector<int> dt (d.cardT, 0);
	for(int i=0; i<d.cardT; ++i){
		//cout << SCIPgetSolVal(scip,sol,xt[i]) << ", ";
		//prod.push_back(static_cast<int>(ceil(SCIPgetSolVal(scip,sol,xt[i]))));
		if(choix==1) dt[i] = static_cast<int>(ceil(d.dt[i]));
		if(choix==2) dt[i] = static_cast<int>(floor(d.dt[i]));
	}
	return dt;

}

float roundd(double var, int nbdec)
{
    // 37.66666 * 100 =3766.66
    // 3766.66 + .5 =3767.16    for rounding off value
    // then type cast to int so value is 3767
    // then divided by 100 so the value converted into 37.67
	float value;
    if(var < 0.0) value = (int)(var * pow(10,nbdec) - 0.5);
	else value = (int)(var * pow(10,nbdec) + 0.5);
    return (float) value / pow(10,nbdec);
}