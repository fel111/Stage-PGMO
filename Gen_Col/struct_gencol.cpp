#include <vector>
#include <iostream>
#include "struct_gencol.h"
using namespace std;


float cl(int k, int l, structGenCol const& sGC){
    return sGC.d.valbpt[0][k];
}

void affK_l(structGenCol const& sGC){
	cout << "--- K_l ---"<<endl;
	for(int l=0; l<sGC.L.size(); ++l){
		cout << l << " : ";
		for(const auto& k : sGC.K_l[l]){
			cout << k << " ";
		}
		cout << endl;
	}
}

void affL_t(structGenCol const& sGC){
	cout << "--- L_t ---"<<endl;
	for(int t=0; t<sGC.d.cardT; ++t){
		cout << t << " : ";
		for(const auto& l : sGC.L_t[t]){
			cout << l << " ";
		}
		cout << endl;
	}
}

void affAllSet(structGenCol const& sGC){
	cout << "--- FeasibleSets ---"<<endl;
	for(int l=0; l<sGC.L.size(); ++l){
		cout << l << " : ";
		for(const auto& t : sGC.L[l].tasksList){
			cout << t << " ";
		}
		cout << endl;
	}
}

// renvoie vrai si l'ensemble l n'est pas déjà présent, faux sinon
bool checkSet(feasibleSet const& l, structGenCol const& sGC){
    for(int i=0; i<sGC.L.size(); ++i){
        if(sGC.L[i].tasksList.size() == l.tasksList.size()){
            int cpt = 0;
            while((sGC.L[i].tasksList[cpt] == l.tasksList[cpt])&&(cpt<l.tasksList.size())) ++cpt;
            if(cpt == l.tasksList.size()) return false;
        }
    }
    return true;
}

void addSetK_l(feasibleSet const& l, structGenCol & sGC){
    vector<int> k;
    for(int p=0; p<(sGC.d.nb_bp[0]/2); ++p){
        if((sGC.d.bpt[0][p*2+1] - l.energyDemand) >= -sGC.p.qmax) k.push_back(p*2);
        if(sGC.d.bpt[0][p*2] - l.energyDemand <= sGC.p.qmax) k.push_back(p*2+1);
    }
    sGC.K_l.push_back(k);
}

void addA_il(feasibleSet const& l, structGenCol & sGC){
    for(int j=0; j<sGC.d.cardJ; ++j){
        int match=0;
        for(const auto& task : l.tasksList){
            if(j==task) match=1;
        }
        sGC.a_il[j].push_back(match);
    }
}

void addL_t(feasibleSet const& l, structGenCol & sGC){
	for(int t=0; t<sGC.d.cardT; ++t){
		if((l.releaseTime<=t)&&(t<l.deadLine)) sGC.L_t[t].push_back(l.id);
	}
}