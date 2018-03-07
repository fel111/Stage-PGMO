#include <iostream>
#include <vector>
#include <chrono>
#include <limits>
#include <deque>

using namespace std;

float infini = numeric_limits<float>::infinity();

float func_fks(vector< vector<float> > f_ks,int k,int s,vector<int> bk){ //renvoie la valeur de f_ks
	if ((s<0) || (s>bk[k-1])) return infini;
	else return f_ks[k][s];
}

float g_ik(int i,int k,vector< vector<float> > f_ks,vector< vector<float> > p_ki,int tau,vector<int> bk){ //renvoie la valeur de g_ik(tau)
	return p_ki[k][i]*(tau)+func_fks(f_ks,k+1,tau,bk);
}

/*float p_t(int t,vector<int> nb_bp,vector<vector <int> > bp_ki,vector<vector <int> > p_ki,vector<vector <int> > f_ki, int x){
	for (int i=0; i<nb_bp[t]; ++i){
		if((x>=bp_ki[t][i])&&(x<bp_ki[t][i+1])) 
	}
}*/

void lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector<vector <int> > bp_ki,vector< vector<float> > p_ki,vector< vector<float> > f_ki,int q_max,int q0){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();

	vector<int> bk (cardT, q_max);
	vector< vector<float> > F_ks (cardT, vector<float> (bk,0.0));
	vector< vector<float> > x_ks (cardT, vector<float> (bk,0.0));

	//calcul F_Ts
	for (int s=0; s<bk; ++s){
		if(s<=dt[cardT-1]) F_ks[cardT-1][s] = p_ki[cardT-1][dt[cardT-1]-s];
		else F_ks[cardT-1][s] = infini;
	}
	
}






int main(){
	cout << "test" << endl;
	return 0;	
}
