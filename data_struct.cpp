#include <vector>
#include "data_struct.h"
using namespace std;


void modifPWL(data &d, vector<int> rt){
	for(int t=0; t<d.cardT; ++t){
		for(int i=0; i<d.nb_bp[t]; ++i){
			d.bpt[t][i] += rt[t];
			//if(rt[t] > 0) d.bpt[t][i] += rt[t];
			//else if(rt[t]<0) d.bpt[t][i] -= rt[t]; 
		}
	}
}