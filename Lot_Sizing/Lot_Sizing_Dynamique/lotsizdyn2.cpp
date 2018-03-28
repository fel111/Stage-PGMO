#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <deque>
#include "lotsizdyn2.h"
#include "data_struct.h"

using namespace std;

float infini_float = numeric_limits<float>::max();
int infini_int = numeric_limits<int>::max();

float func_fks(vector< vector<float> > f_ks,int k,int s,int bk){ //renvoie la valeur de f_ks
	//cout << "func FKS : "<<k<< " ," << s<<endl;	
	if ((s<0) || (s>bk)) return infini_float;
	else return f_ks[k][s];
}

float p_t(int t,vector<int> nb_bp,vector<vector <int> > bpt,vector< vector<float> > valbpt,vector<vector <float> > pente, int x){
	//cout << "x: "<<x<<endl;	
	if (x<bpt[t][0]) return 0.0;	
	for (int i=0; i<nb_bp[t]; ++i){
		//cout<<bpt[t][i]<<"<"<<x<<"<="<<bpt[t][i+1]<<endl;
		if((x>bpt[t][i])&&(x<=bpt[t][i+1])){
			//cout << bpt[t][i] << endl;
			//cout << "valreel" <<  pente[t][i]*x+valbpt[t][i]<<endl;
			return pente[t][i]*(x-bpt[t][i])+valbpt[t][i];
			//return pente[t][i]*x+valbpt[t][i];
		}	
	}
	return infini_float;
}

float g_ki(int i,int k,vector< vector<float> > f_ks,vector<vector <float> > pente, int tau,int bk){
	return pente[k][i]*tau+func_fks(f_ks,k+1,tau,bk);
	
}

//a inclure dans un header ensuite
void lecture_pwd(string file, vector<vector<int> >& bpt, vector<vector<float> >& valbpt, vector<vector<float> >& pente){

	//char temp;
	//vector< vector<int> > conflitstemp;
	
	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		int bp;
		float valbp;
		float p;
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> bp >> valbp >> p;
			bpt[t].push_back(bp);
			valbpt[t].push_back(valbp);
			pente[t].push_back(p);
			//cout << "OK" << endl;
		}
		fichier.close();
		//cout << " done "<<endl;
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}

}

//a inclure dans un header ensuite
void lecture_demande(string file, vector<int> & dt){

	//char temp;
	//vector< vector<int> > conflitstemp;
	
	ifstream fichier(file, ios::in); //ouverture du fichier
	if(fichier){  // si l'ouverture fonctionne
		string ligne;
		int t;
		int d;
		while(getline(fichier,ligne)){ // tant que l'on peut lire une ligne
			istringstream iss (ligne);
			iss >> t >> d;
			dt.push_back(d);
			//cout << "dt["<<t<<"] : "<<dt[t]<<endl;
		}
		fichier.close();
	}
	else{
		cout << "erreur lecture fichier" << endl;
	}
}
/*
void lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector< vector<float> > valbpt,vector< vector<int> > bpt,vector<vector <float> > pente,int q_max){//vector< vector<float> > f_ki,vector< vector<float> > g_ki){//vector<int> q0){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();
		
	vector< vector<float> > f_ki (cardT, vector<float> (nb_bp[0],0.0));
	for(int k=0; k<cardT; ++k){
		for(int i=0; i<nb_bp[k]; ++i){
			f_ki[k][i] = valbpt[k][i] - pente[k][i]*bpt[k][i];
			cout << f_ki[k][i] <<endl;
		}
	}

	//vector<int> bk (cardT, q_max);
	vector< vector<float> > F_ks (cardT, vector<float> (q_max+1,0.0));
	vector< vector<int> > x_ks (cardT, vector<int> (q_max+1,0.0));

	//calcul F_Ts
	//cout <<"F["<<cardT-1<<"] = [";
	for (int s=0; s<=q_max; ++s){
		if(s<=dt[cardT-1]) F_ks[cardT-1][s] = p_t(cardT-1,nb_bp,bpt,valbpt,pente,dt[cardT-1]-s);
		else F_ks[cardT-1][s] = infini_float;
		//cout << F_ks[cardT-1][s] <<" ";
	}
	//cout << "]"<<endl;
	
	int k = cardT-1;
	//deque<vector< vector<float> > > G_ki;
	//deque<vector<vector<int> > > x_kis;
	while (k!=0){
		--k;
		//cout << "k : "<<k<<endl;
		vector< vector<float> >Gk_i;
		vector< vector<int> >xk_is;
		for(int i=0; i<nb_bp[k]; ++i){
			//cout << "i : "<<i<<endl;
			vector<float> Gki ;
			vector<int> xki_s ;
			int s=0;
			bool seq = false;
			while((s<=q_max)&&(seq!=true)){
				//cout <<"while s : " << s << endl;
				int lb = max(bpt[k][i]+1+s-dt[k],0);
				//int lb = max(bpt[k][i]+s-dt[k],0);
				//int ub = min(bpt[k][i+1]-1+s-dt[k],q_max);
				int ub = min(bpt[k][i+1]+s-dt[k],q_max);
				//cout << "lb, ub : "<<lb <<" "<<ub<<endl;
				if(lb>ub){
					//cout << "SEQUENCE IMPOSSIBLE NUMERO "<<s<<endl;
					//cout <<lb<<" "<<ub<<endl;
					Gki.push_back(infini_float);
					++s;					
					xki_s.push_back(infini_int);
					//xks[i] = infini;
				}
				else{
					//cout << "PREMIER S POSSIBLE NUMERO "<<s<<endl;
					//cout <<lb<<" "<<ub<<endl;
					
					seq = true; // on a un premier s possible
					deque<int> u = {ub};
					--ub;
					while(lb<=ub){
						//cout << "while(lb<=ub)"<< endl;
						if(g_ki(i,k,F_ks,pente,ub,q_max)<g_ki(i,k,F_ks,pente,u[0],q_max)){
							u.push_front(ub);
							--ub;
						}
						else{
							--ub;
						}
					} // ici, les u sont triés selon g(u1) < g(u2) < ...
											
					Gki.push_back(g_ki(i,k,F_ks,pente,u[0],q_max));
					if(dt[k]-s+u[0]>=0) xki_s.push_back(dt[k]-s+u[0]);
					else xki_s.push_back(infini_int);
					++s;
					while(s<=q_max){ // on parcours les s et on met à jour la sequence u
							
						//if((u.size()!=0)&&((bpt[k][i]+1+s-dt[k])==u[0])){
						if((u.size()!=0)&&((bpt[k][i]+s-dt[k])==u[0])){
							u.pop_front();
						}
						if((bpt[k][i+1]+s-dt[k])<=q_max){
							while((u.size()!=0)&&(g_ki(i,k,F_ks,pente,u.back(),q_max)>=g_ki(i,k,F_ks,pente,bpt[k][i+1]+s-dt[k],q_max))){
								u.pop_back();
							}
							u.push_back(bpt[k][i+1]+s-dt[k]);
						}
						//cout << "SEQUENCE U MAJ NUMERO "<<s<<endl;
						//for(int j=0;j<u.size();++j){ cout << u[j] << " ";}
						//cout << endl;
						if(u.size()!=0){
							Gki.push_back(g_ki(i,k,F_ks,pente,u[0],q_max));
							if(dt[k]-s+u[0]>=0) xki_s.push_back(dt[k]-s+u[0]);
							else xki_s.push_back(infini_int);
						}
						else{
							Gki.push_back(infini_float);
							xki_s.push_back(infini_int);
						}
						++s;
					}
				}
			}
			
			Gk_i.push_back(Gki);
			//cout << "xki_s pour k et i fixé =" <<endl;
			//for(int j=0;j<xki_s.size();++j){ cout << xki_s[j] << " ";}
			//cout << endl;
			xk_is.push_back(xki_s);
		}
		//G_ki.push_front(Gk_i);
		//x_kis.push_front(xk_is);
		//cout <<"x["<<k<<"] = [";
		for(int s=0;s<=q_max;++s){
			float fGauche;
			fGauche = func_fks(F_ks,k+1,s-dt[k],q_max);
			float minFdroite = f_ki[k][0]+pente[k][0]*(dt[k]-s)+Gk_i[0][s];
			int imin = 0;			
			for(int i=1;i<nb_bp[k];++i){
				if((f_ki[k][i]+pente[k][i]*(dt[k]-s)+Gk_i[i][s]>=0)&&(f_ki[k][i]+pente[k][i]*(dt[k]-s)+Gk_i[i][s]<minFdroite)){
					minFdroite = f_ki[k][i]+pente[k][i]*(dt[k]-s)+Gk_i[i][s];
					imin = i;
				}			
			}
			
			if(fGauche<minFdroite){
				F_ks[k][s] = fGauche;
				x_ks[k][s] = 0;
				//cout <<"F["<<k<<"]["<<s<<"] = "<<F_ks[k][s] <<endl;
				//cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
			}
			else{
				F_ks[k][s] = minFdroite;
				//cout << "imin : "<<imin<<endl;
				//cout << "x_kis : "<<x_kis[0][imin][s]<<endl;
				x_ks[k][s] = xk_is[imin][s];	
				//x_ks[k][s] = x_kis[0][imin][s];
				//cout <<"F["<<k<<"]["<<s<<"] = "<<F_ks[k][s] <<endl;
				//cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
			}
			//cout << x_ks[k][s] << " ";			

		}
		//cout << "]"<<endl;
		
	}
	for(int k=0; k<cardT; ++k){
		for(int s=0; s<=q_max; ++s){
			cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
		}
	}
	
	cout << "F[0][0] = " <<F_ks[0][0] <<endl;
	//'''determination de xt'''    CALCUL DES XT ET DU STOCK BATTERIE
	vector<int> stock (cardT,0);// = q0;
	vector<int> xt (cardT,0);
	xt[0] = x_ks[0][0];
	stock[0] = xt[0] - dt[0];
	xt[1] = x_ks[1][stock[0]];
	for (int z=1; z<cardT-2; ++z){		
		stock[z]=(stock[z-1]+xt[z]-dt[z]);
		xt[z+1] = (x_ks[z+1][stock[z]]);
		cout << "stock" << z << " : " << stock[z] << endl;
		cout << "xt" << z << " : " << xt[z+1] << endl;
		
	}
	stock[cardT-2] = (stock[cardT-3]+xt[cardT-2]-dt[cardT-2]);
	xt[cardT-1] = dt[cardT-1] - stock[cardT-2];
	stock[cardT-1] = (stock[cardT-2]+xt[cardT-1]-dt[cardT-1]);

	for(int i=0;i<cardT;++i){
		cout << "stock" << i << " : " << stock[i] << endl;
	}	
	for(int i=0;i<cardT;++i){
		cout << "xt" << i << " : " << xt[i] << endl;	
	}
	float cout_tot = 0;
	for(int i=0;i<cardT;++i){
		cout_tot += cout_reel(i,nb_bp,bpt,valbpt,pente,xt[i]);
	}
	//cout << "cout total : "<<cout_tot<<endl;
	cout << "cout total : "<<F_ks[0][0]<<endl;
	//return ;
}*/



vector<int> lotsizdyn(data d){//vector< vector<float> > f_ki,vector< vector<float> > g_ki){//vector<int> q0){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();
		
	vector< vector<float> > f_ki (d.cardT, vector<float> (d.nb_bp[0],0.0));
	for(int k=0; k<d.cardT; ++k){
		for(int i=0; i<d.nb_bp[k]; ++i){
			f_ki[k][i] = d.valbpt[k][i] - d.pente[k][i]*d.bpt[k][i];
			//cout << f_ki[k][i] <<endl;
		}
	}

	//vector<int> bk (d.cardT, d.Q);
	vector< vector<float> > F_ks (d.cardT, vector<float> (d.Q+1,0.0));
	vector< vector<int> > x_ks (d.cardT, vector<int> (d.Q+1,0.0));

	//calcul F_Ts
	//cout <<"F["<<d.cardT-1<<"] = [";
	for (int s=0; s<=d.Q; ++s){
		if(s<=d.dt[d.cardT-1]) F_ks[d.cardT-1][s] = p_t(d.cardT-1,d.nb_bp,d.bpt,d.valbpt,d.pente,d.dt[d.cardT-1]-s);
		else F_ks[d.cardT-1][s] = infini_float;
		//cout << F_ks[d.cardT-1][s] <<" ";
	}
	//cout << "]"<<endl;
	
	int k = d.cardT-1;
	//deque<vector< vector<float> > > G_ki;
	//deque<vector<vector<int> > > x_kis;
	while (k!=0){
		--k;
		//cout << "k : "<<k<<endl;
		vector< vector<float> >Gk_i;
		vector< vector<int> >xk_is;
		for(int i=0; i<d.nb_bp[k]; ++i){
			//cout << "i : "<<i<<endl;
			vector<float> Gki ;
			vector<int> xki_s ;
			int s=0;
			bool seq = false;
			while((s<=d.Q)&&(seq!=true)){
				//cout <<"while s : " << s << endl;
				int lb = max(d.bpt[k][i]+1+s-d.dt[k],0);
				//int lb = max(d.bpt[k][i]+s-d.dt[k],0);
				//int ub = min(d.bpt[k][i+1]-1+s-d.dt[k],d.Q);
				int ub = min(d.bpt[k][i+1]+s-d.dt[k],d.Q);
				//cout << "lb, ub : "<<lb <<" "<<ub<<endl;
				if(lb>ub){
					//cout << "SEQUENCE IMPOSSIBLE NUMERO "<<s<<endl;
					//cout <<lb<<" "<<ub<<endl;
					Gki.push_back(infini_float);
					++s;					
					xki_s.push_back(infini_int);
					//xks[i] = infini;
				}
				else{
					//cout << "PREMIER S POSSIBLE NUMERO "<<s<<endl;
					//cout <<lb<<" "<<ub<<endl;
					
					seq = true; // on a un premier s possible
					deque<int> u = {ub};
					--ub;
					while(lb<=ub){
						//cout << "while(lb<=ub)"<< endl;
						if(g_ki(i,k,F_ks,d.pente,ub,d.Q)<g_ki(i,k,F_ks,d.pente,u[0],d.Q)){
							u.push_front(ub);
							--ub;
						}
						else{
							--ub;
						}
					} // ici, les u sont triés selon g(u1) < g(u2) < ...
											
					Gki.push_back(g_ki(i,k,F_ks,d.pente,u[0],d.Q));
					if(d.dt[k]-s+u[0]>=0) xki_s.push_back(d.dt[k]-s+u[0]);
					else xki_s.push_back(infini_int);
					++s;
					while(s<=d.Q){ // on parcours les s et on met à jour la sequence u
							
						//if((u.size()!=0)&&((d.bpt[k][i]+1+s-d.dt[k])==u[0])){
						if((u.size()!=0)&&((d.bpt[k][i]+s-d.dt[k])==u[0])){
							u.pop_front();
						}
						if((d.bpt[k][i+1]+s-d.dt[k])<=d.Q){
							while((u.size()!=0)&&(g_ki(i,k,F_ks,d.pente,u.back(),d.Q)>=g_ki(i,k,F_ks,d.pente,d.bpt[k][i+1]+s-d.dt[k],d.Q))){
								u.pop_back();
							}
							u.push_back(d.bpt[k][i+1]+s-d.dt[k]);
						}
						//cout << "SEQUENCE U MAJ NUMERO "<<s<<endl;
						//for(int j=0;j<u.size();++j){ cout << u[j] << " ";}
						//cout << endl;
						if(u.size()!=0){
							Gki.push_back(g_ki(i,k,F_ks,d.pente,u[0],d.Q));
							if(d.dt[k]-s+u[0]>=0) xki_s.push_back(d.dt[k]-s+u[0]);
							else xki_s.push_back(infini_int);
						}
						else{
							Gki.push_back(infini_float);
							xki_s.push_back(infini_int);
						}
						++s;
					}
				}
			}
			
			Gk_i.push_back(Gki);
			//cout << "xki_s pour k et i fixé =" <<endl;
			//for(int j=0;j<xki_s.size();++j){ cout << xki_s[j] << " ";}
			//cout << endl;
			xk_is.push_back(xki_s);
		}
		//G_ki.push_front(Gk_i);
		//x_kis.push_front(xk_is);
		//cout <<"x["<<k<<"] = [";
		for(int s=0;s<=d.Q;++s){
			float fGauche;
			fGauche = func_fks(F_ks,k+1,s-d.dt[k],d.Q);
			float minFdroite = f_ki[k][0]+d.pente[k][0]*(d.dt[k]-s)+Gk_i[0][s];
			int imin = 0;			
			for(int i=1;i<d.nb_bp[k];++i){
				if((f_ki[k][i]+d.pente[k][i]*(d.dt[k]-s)+Gk_i[i][s]>=0)&&(f_ki[k][i]+d.pente[k][i]*(d.dt[k]-s)+Gk_i[i][s]<minFdroite)){
					minFdroite = f_ki[k][i]+d.pente[k][i]*(d.dt[k]-s)+Gk_i[i][s];
					imin = i;
				}			
			}
			
			if(fGauche<minFdroite){
				F_ks[k][s] = fGauche;
				x_ks[k][s] = 0;
				//cout <<"F["<<k<<"]["<<s<<"] = "<<F_ks[k][s] <<endl;
				//cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
			}
			else{
				F_ks[k][s] = minFdroite;
				//cout << "imin : "<<imin<<endl;
				//cout << "x_kis : "<<x_kis[0][imin][s]<<endl;
				x_ks[k][s] = xk_is[imin][s];	
				//x_ks[k][s] = x_kis[0][imin][s];
				//cout <<"F["<<k<<"]["<<s<<"] = "<<F_ks[k][s] <<endl;
				//cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
			}
			//cout << x_ks[k][s] << " ";			

		}
		//cout << "]"<<endl;
		
	}
	/*for(int k=0; k<d.cardT; ++k){
		for(int s=0; s<=d.Q; ++s){
			cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
		}
	}*/
	
	//cout << "F[0][0] = " <<F_ks[0][0] <<endl;
	//'''determination de xt'''    CALCUL DES XT ET DU STOCK BATTERIE
	vector<int> rt (d.cardT,0);
	vector<int> stock (d.cardT,0);// = q0;
	vector<int> xt (d.cardT,0);
	xt[0] = x_ks[0][0];
	stock[0] = xt[0] - d.dt[0];
	//cout << "stock" << stock[0] << endl;
	xt[1] = x_ks[1][stock[0]];
	for (int z=1; z<d.cardT-2; ++z){		
		stock[z]=(stock[z-1]+xt[z]-d.dt[z]);
		xt[z+1] = (x_ks[z+1][stock[z]]);
		//cout << "stock" << z << " : " << stock[z] << endl;
		//cout << "xt" << z << " : " << xt[z+1] << endl;
		
	}
	stock[d.cardT-2] = (stock[d.cardT-3]+xt[d.cardT-2]-d.dt[d.cardT-2]);
	xt[d.cardT-1] = d.dt[d.cardT-1] - stock[d.cardT-2];
	stock[d.cardT-1] = (stock[d.cardT-2]+xt[d.cardT-1]-d.dt[d.cardT-1]);

	rt[0] = stock[0] - d.s0;
	/*cout << "rt" << endl;
	cout << rt[0] << endl;*/
	for(int t=1; t<d.cardT; ++t){
		rt[t] = stock[t] - stock[t-1];
		//cout << rt[t] << endl;
	}

	/*for(int i=0;i<d.cardT;++i){
		cout << "stock" << i << " : " << stock[i] << endl;
	}	
	for(int i=0;i<d.cardT;++i){
		cout << "xt" << i << " : " << xt[i] << endl;	
	}*/
	/*float cout_tot = 0;
	for(int i=0;i<d.cardT;++i){
		cout_tot += cout_reel(i,d.nb_bp,d.bpt,d.vald.bpt,d.pente,xt[i]);
	}*/
	//cout << "cout total : "<<cout_tot<<endl;
	cout << "cout total : "<<F_ks[0][0]<<endl;
	return rt;
}

/*
int main(){
	//creation donnees test
	int cardT = 10; 
	int qmax = 20;
	vector<int> d.dt = {5,10,15,20,25,12,45,24,12,65};
	vector<vector<int> > d.bpt;
	vector<vector<float> > d.vald.bpt;
	vector<vector<float> > d.pente;
	vector<int> nb_bp = {3,3,3,3,3,3,3,3,3,3};
	for(int i=0;i<cardT;++i){
		vector<float> temp_pente = {3.0,2.0,1.0};
		vector<int> temp_d.bpt = {0,10,20,infini_int};
		vector<float> temp_d.valbpt = {0.0,30.0,50.0};
		pente.push_back(temp_pente);
		bpt.push_back(temp_bpt);
		valbpt.push_back(temp_valbpt);
	}
	int cardT = 3; 
	int qmax = 5;
	vector<int> dt = {5,10,10};
	vector<vector<int> > bpt;
	vector<vector<float> > valbpt;
	vector<vector<float> > pente;
	vector<int> nb_bp = {2,2,2};
	for(int i=0;i<cardT;++i){
		vector<float> temp_pente = {2.0,1.0};
		vector<int> temp_bpt = {0,10,infini_int};
		vector<float> temp_valbpt = {0.0,20.0};
		pente.push_back(temp_pente);
		bpt.push_back(temp_bpt);
		valbpt.push_back(temp_valbpt);
	}
	
	int cardT = 100;
	int qmax = 50;
	vector<int> nb_bp (cardT, 10);
	vector<vector<int> > bpt;
	vector<vector<float> > valbpt;
	vector<vector<float> > pente;
	for(int i=0; i<cardT; ++i){
		vector<int> bp;
		vector<float> valbp;
		vector<float> p;
		bpt.push_back(bp);
		valbpt.push_back(valbp);
		pente.push_back(p);	
	}
	vector<int> dt;
	
	lecture_demande("demande_100_200", dt);
	lecture_pwd("pwd_100_10_3_100", bpt, valbpt, pente);

	for(int k=0; k<cardT;++k){
		bpt[k].push_back(infini_int);
	}
	for(int k=0; k<cardT;++k){
		for(int i=0; i<nb_bp[k]; ++i){
			bpt[][]
			cout << "bpt["<<k<<"]["<<i<<"] = "<<bpt[k][i]<<endl;
			cout << "valbpt["<<k<<"]["<<i<<"] = "<<valbpt[k][i]<<endl;
			cout << "pente["<<k<<"]["<<i<<"] = "<<pente[k][i]<<endl;
		}
	}

	chrono::steady_clock::time_point begin = chrono::steady_clock::now();
	lotsizing(cardT,dt,nb_bp,valbpt,bpt,pente,qmax);
	chrono::steady_clock::time_point end= chrono::steady_clock::now();
	cout << "Temps = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() <<endl;
	cout << "Temps = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() /1000000.0 <<endl;
	for(int i=0; i<xt.size(); ++i){
		cout << "xt"<<i<<" : "<<xt[i] << endl;
	}
	return 0;	
}*/