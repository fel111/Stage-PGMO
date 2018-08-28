#include <vector>
#include <iostream>
#include <limits>
#include "struct_gencol.h"
using namespace std;
float infinif = numeric_limits<float>::max();

float cl(int k, int l, structGenCol const& sGC){
    return sGC.d->valbpt[0][k];
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
	for(int t=0; t<sGC.d->cardT; ++t){
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

// renvoie l'id de l'ensemble si deja present, -1 sinon
int checkSet(feasibleSet const& l, structGenCol const& sGC){
    for(int i=0; i<sGC.L.size(); ++i){
        if(sGC.L[i].tasksList.size() == l.tasksList.size()){
            int cpt = 0;
			/*cout << "1:"<<sGC.L[i].tasksList[cpt] << endl;
			cout << "2:"<<l.tasksList[cpt] << endl;
			cout << "3:"<<l.tasksList.size() << endl;
			bool suite = false;
			if((sGC.L[i].tasksList[cpt] == l.tasksList[cpt])){
				if(cpt<l.tasksList.size()){
					suite = true;
				}
			}*/
			bool stop = true;
			if(cpt<l.tasksList.size()){
				if((sGC.L[i].tasksList[cpt] == l.tasksList[cpt])){
					stop = false;
				}
			}
			while(!stop){
				++cpt;
				if(cpt>=l.tasksList.size()){
					stop = true;
				}
				else if((sGC.L[i].tasksList[cpt] != l.tasksList[cpt])){
						stop = true;
				}
			}
			
            if(cpt == l.tasksList.size()) return i;
        }
    }
    return -1;
}

/*void addSetK_l(feasibleSet const& l, structGenCol & sGC){
    vector<int> k;
    for(int p=0; p<(sGC.d->nb_bp[0]/2); ++p){
        if((sGC.d->bpt[0][p*2+1] - l.energyDemand) >= -sGC.d->Q) k.push_back(p*2);
        if(sGC.d->bpt[0][p*2] - l.energyDemand <= sGC.d->Q) k.push_back(p*2+1);
    }
    sGC.K_l.push_back(k);
}*/

// correction 
void addSetK_l(feasibleSet const& l, structGenCol & sGC){
    vector<int> k;
	vector<int> pl;
    for(int p=0; p<(sGC.d->nb_bp[0]/2); ++p){
		if((sGC.d->bpt[0][p*2] - l.energyDemand <= sGC.d->Q)&&(l.energyDemand - sGC.d->bpt[0][p*2+1] <= sGC.d->Q)){
			k.push_back(p*2);
			k.push_back(p*2+1);
			pl.push_back(1);
		}
		else pl.push_back(0);
    }
    sGC.K_l.push_back(k);
	sGC.P_l.push_back(pl);
}


void addA_il(feasibleSet const& l, structGenCol & sGC){
    for(int j=0; j<sGC.d->cardJ; ++j){
        int match=0;
        for(const auto& task : l.tasksList){
            if(j==task) match=1;
        }
        sGC.a_il[j].push_back(match);
    }
}

/*void addA_il(feasibleSet const& l, structGenCol & sGC){
	for(const auto& task : l.tasksList){
		sGC.a_il[task].push_back(l.id);
	}
}*/

void addL_t(feasibleSet const& l, structGenCol & sGC){
	for(int t=0; t<sGC.d->cardT; ++t){
		if((l.releaseTime<=t)&&(t<l.deadLine)) sGC.L_t[t].push_back(l.id);
	}
}

void addP_0(structGenCol & sGC){
	vector<int> p0;
	for(int p=0; p<(sGC.d->nb_bp[0]/2); ++p){
		if(sGC.d->bpt[0][p*2] <= sGC.d->Q) p0.push_back(1);
		else p0.push_back(0);
	}
	sGC.P_0 = p0;
}

void modifPwlCplex(structGenCol & sGC){ // a rendre generique selon les pwl
	vector<int> nbtemp;
	for(int i=0;i<sGC.d->cardT;++i){
		sGC.d->bpt[i] = {0.0,5.0,100.0};
		sGC.d->valbpt[i] = {0.0,10.0,152.5};
		sGC.d->pente[i] = {2.0,1.5};
		nbtemp.push_back(2);
	}
	sGC.d->nb_bp = nbtemp;
}

float p_t(int t, float x, structGenCol const& sGC){
	//cout << "x: "<<x<<endl;	
	if (x<sGC.d->bpt[t][0]) return 0.0;	
	for (int i=0; i<sGC.d->nb_bp[t]; ++i){
		//cout<<bpt[t][i]<<"<"<<x<<"<="<<bpt[t][i+1]<<endl;
		if((sGC.d->bpt[t][i]<=x)&&(x<sGC.d->bpt[t][i+1])){
			//cout << bpt[t][i] << endl;
			//cout << "valreel" <<  pente[t][i]*x+valbpt[t][i]<<endl;
			return sGC.d->pente[t][i]*(x-sGC.d->bpt[t][i])+sGC.d->valbpt[t][i];
			//return pente[t][i]*x+valbpt[t][i];
		}	
	}
	return infinif;
}