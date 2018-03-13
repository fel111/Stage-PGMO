#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <deque>

using namespace std;

float infini = numeric_limits<float>::infinity();

/*float func_fks(vector< vector<float> > f_ks,int k,int s,vector<int> bk){ //renvoie la valeur de f_ks
	if ((s<0) || (s>bk[k-1])) return infini;
	else return f_ks[k][s];
}*/

float p_t(int t,vector<int> nb_bp,vector<vector <int> > bpt,vector<vector <int> > valbpt,vector<vector <float> > pente, float x){
	if (x==0.0) return 0.0;	
	for (int i=0; i<nb_bp[t]; ++i){
		if((x>bpt[t][i])&&(x<=bpt[t][i+1])) 
			return pente[t][i]*x+valbpt[t][i];
	}
}

/*float g_ki(int k, int i, int tau,){
	return pente[k][i]*tau+func_fks(k+1,tau)
	
}*/

/*float G_ik(int i, int k, int eta, vector< vector<float> > bpt,int & u0, bool upd,vector<int> dt,int s,vector<vector <int> > valbpt,int q_max){
	float lb = max(bpt[k][i]+1+s-dt[k],0);
	float up = min(bpt[k][i+1]+s-dt[k],q_max);
	if(lb>ub) return infini;
	float min = g_ki[k][lb];
	int indmin = 0;
	for(int j=lb+1; j<lb+up; ++j){
		if(min>g_ki[k][j]){ min = g_ki[k][j]; indmin=j;}	
	}
	if(upd) u0 = indmin;
	return min;

}*/
/*float p_t(int t,vector<int> nb_bp,vector<vector <int> > bbpt,vector<vector <int> > bpt,vector<vector <int> > f_ki, int x){
	for (int i=0; i<nb_bp[t]; ++i){
		if((x>=bbpt[t][i])&&(x<bbpt[t][i+1])) 
	}
}*/

vector<float> lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector<vector <int> > valbpt,vector< vector<float> > bpt,int q_max){//vector< vector<float> > f_ki,vector< vector<float> > g_ki){//vector<int> q0){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();

	vector<int> bk (cardT, q_max);
	vector< vector<float> > F_ks (cardT, vector<float> (q_max,0.0));
	vector< vector<float> > x_ks (cardT, vector<float> (q_max,0.0));

	//calcul F_Ts
	for (int s=0; s<bk[cardT-1]; ++s){
		if(s<=dt[cardT-1]) F_ks[cardT-1][s] = p_t(cardT-1,nb_bp,bpt,valbpt,pente,dt[cardT-1]-s);
		else F_ks[cardT-1][s] = infini;
	}
	
	int k = cardT-2;
	while (k!=0){
		//vector<float> G_ik (nb_bp[k],0.0);
		/*for(int s=0; s<bk[cardT-2]; ++s){
			int u0;
			float minval = valbpt[k][0] + bpt[k][0]*(dt[k]-s) + G_ik(0,k,(s-dt[k]),g_ki,u0,true,dt,s,valbpt,q_max);
			for(int i=1; i<nb_bp[k]; ++i){
				if (minval > (f_ki[k][i] + bpt[k][i]*(dt[k]-s)+G_ik(i,k,(s-dt[k]),g_ki,u0,false,dt,s,valbpt,q_max))) minval = (f_ki[k][i] + bpt[k][i]*(dt[k]-s)+G_ik(i,k,(s-dt[k]),g_ki,u0,true,dt,s,valbpt,q_max));	
			}
			F_ks[k][s] = min(F_ks[k+1][s-dt[k]],minval);
			x_ks[k][s] = dt[k]-s+u0;
		}*/
		
		for(int s=0; s<bk[cardT-2]; ++s){
			vector<float> Gki (nb_bp[k],infini);
			vector<float> xks (nb_bp[k],infini);
			for(int i=0; i<nb_bp[k]-1; ++i){
				float lb = max(bpt[k][i]+1+s-dt[k],0);
				float up = min(bpt[k][i+1]+s-dt[k],q_max);
				if(lb>ub){
					Gki[i] = infini;
					xks[i] = infini;
				}
				else{
					Gki[i] = p_t(k,nb_bp,bpt,valbpt,pente,lb);
				}
			}
		}
		--k;
	}
	//'''determination de xt'''    CALCUL DES XT ET DU STOCK BATTERIE
	vector<int> stock (cardT,0);// = q0;
	vector<float> xt (cardT,0.0);
	xt[0] = x_ks[0][0];
	for (int z=0; z<cardT-1; ++z){		
		stock[z]=(xt[z]-dt[z]);
		cout << "stock" << z << " : " << stock[z] << endl;
		xt[z+1] = (x_ks[z+1][stock[z]]);
	}
	stock[cardT-1] = (xt[cardT-1]-dt[cardT-1]);
	return xt;
}






int main(){
	//creation donnees test
	/*int cardT = 10; 
	int q_max = 100;
	vector<int> bp = {0,20,40,60,80};
	vector<float> p={1.0,0.6,0.8,0.7,0.6};
	vector<int> dt = {5,10,15,20,25,30,35,40,45,50};
	vector<vector<int> > valbpt;
	vector<int> nb_bp;
	vector<vector<float> > bpt;
	vector<vector<float> > f_ki;
	vector<vector<float> > g_ki;
	vector<float> f = {0.0};
	
	for(int i=0; i<cardT; ++i){ 
		valbpt.push_back(bp);
		bpt.push_back(p);
		nb_bp.push_back(valbpt[i].size());
	}
	for(int i=0; i<nb_bp[i]-1; ++i){
		float temp = p[i]*bp[i+1]+f[i]-p[i+1]*bp[i+1];
		f.push_back(temp);
	}
	for(int i=0; i<cardT; ++i){
		f_ki.push_back(f);
		//f_ki[i].push_back(infini);
	}
	g_ki = f_ki;*/
	
	int cardT = 5; 
	int qmax = 20;
	vector<int> dt = {5,10,15,20,25};
	vector<vector<int> > bpt;
	vector<vector<int> > valbpt;
	vector<vector<float> > pente;
	vector<int> nb_bp = {3,3,3,3,3};
	for(int i=0;i<cardT;++i){
		vector<float> temp_pente = {1.0,2.0,1.0};
		vector<int> temp_bpt = {0,10,20};
		vector<int> temp_valbpt = {0,10,30};
		pente.push_back(temp_pente);
		bpt.push_back(temp_bpt);
		valbpt.push_back(temp_valbpt);
	}
	
	vector<float> xt = lotsizing(cardT,dt,nb_bp,val_bpt,bpt,qmax);
	for(int i=0; i<xt.size(); ++i){
		cout << xt[i] << endl;
	}
	return 0;	
}
