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

/*float g_ik(int i,int k,vector< vector<float> > f_ks,vector< vector<float> > p_ki,int tau,vector<int> bk){ //renvoie la valeur de g_ik(tau)
	return p_ki[k][i]*(tau)+func_fks(f_ks,k+1,tau,bk);
}*/

float G_ik(int i, int k, int eta, vector< vector<float> > g_ki,int & u0, bool upd,vector<int> dt,int s,vector<vector <int> > P_ki,int q_max){
	float lb = max(P_ki[k][i]+1+s-dt[k],0);
	float up = min(P_ki[k][i+1]+s-dt[k],q_max);
	float min = g_ki[k][lb];
	int indmin = 0;
	for(int j=lb+1; j<lb+up; ++j){
		if(min>g_ki[k][j]){ min = g_ki[k][j]; indmin=j;}	
	}
	if(upd) u0 = indmin;
	return min;

}
/*float p_t(int t,vector<int> nb_bp,vector<vector <int> > bp_ki,vector<vector <int> > p_ki,vector<vector <int> > f_ki, int x){
	for (int i=0; i<nb_bp[t]; ++i){
		if((x>=bp_ki[t][i])&&(x<bp_ki[t][i+1])) 
	}
}*/

vector<float> lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector<vector <int> > P_ki,vector< vector<float> > p_ki,vector< vector<float> > f_ki,vector< vector<float> > g_ki,int q_max){//vector<int> q0){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();

	vector<int> bk (cardT, q_max);
	vector< vector<float> > F_ks (cardT, vector<float> (q_max,0.0));
	vector< vector<float> > x_ks (cardT, vector<float> (q_max,0.0));

	//calcul F_Ts
	for (int s=0; s<bk[cardT-1]; ++s){
		if(s<=dt[cardT-1]) F_ks[cardT-1][s] = p_ki[cardT-1][dt[cardT-1]-s];
		else F_ks[cardT-1][s] = infini;
	}
	
	int k = cardT-2;
	while (k!=0){
		//vector<float> G_ik (nb_bp[k],0.0);
		for(int s=0; s<bk[cardT-2]; ++s){
			int u0;
			float minval = f_ki[k][0] + p_ki[k][0]*(dt[k]-s)+G_ik(0,k,(s-dt[k]),g_ki,u0,true,dt,s,P_ki,q_max);
			for(int i=1; i<nb_bp[k]; ++i){
				if (minval > (f_ki[k][i] + p_ki[k][i]*(dt[k]-s)+G_ik(i,k,(s-dt[k]),g_ki,u0,false,dt,s,P_ki,q_max))) minval = (f_ki[k][i] + p_ki[k][i]*(dt[k]-s)+G_ik(i,k,(s-dt[k]),g_ki,u0,true,dt,s,P_ki,q_max));	
			}
			F_ks[k][s] = min(F_ks[k+1][s-dt[k]],minval);
			x_ks[k][s] = dt[k]-s+u0;
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
	int cardT = 10; 
	int q_max = 100;
	vector<int> bp = {0,20,40,60,80};
	vector<float> p={1.0,0.6,0.8,0.7,0.6};
	vector<int> dt = {5,10,15,20,25,30,35,40,45,50};
	vector<vector<int> > P_ki;
	vector<int> nb_bp;
	vector<vector<float> > p_ki;
	vector<vector<float> > f_ki;
	vector<vector<float> > g_ki;
	vector<float> f = {0.0};
	
	for(int i=0; i<cardT; ++i){ 
		P_ki.push_back(bp);
		p_ki.push_back(p);
		nb_bp.push_back(P_ki[i].size());
	}
	for(int i=0; i<nb_bp[i]-1; ++i){
		float temp = p[i]*bp[i+1]+f[i]-p[i+1]*bp[i+1];
		f.push_back(temp);
	}
	for(int i=0; i<cardT; ++i){
		f_ki.push_back(f);
		//f_ki[i].push_back(infini);
	}
	g_ki = f_ki;
	
	vector<float> xt = lotsizing(cardT,dt,nb_bp,P_ki,p_ki,f_ki,g_ki,q_max);
	for(int i=0; i<xt.size(); ++i){
		cout << xt[i] << endl;
	}
	return 0;	
}
