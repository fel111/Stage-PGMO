#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <deque>
#include "lotsizdyn2.h"
#include "struct.h"

using namespace std;

constexpr float infini_float = numeric_limits<float>::infinity();
constexpr int infini_int = numeric_limits<int>::max();

float func_fks(vector< vector<float> > const& f_ks,int k,int s,int bk){ //renvoie la valeur de f_ks
	//cout << "func FKS : "<<k<< " ," << s<<endl;	
	if ((s<0) || (s>bk)) return infini_float;
	else return f_ks[k][s];
}

float p_t(int t,vector<int> const& nb_bp,vector<vector <float> > const& bpt, vector< vector<float> > const& valbpt,vector<vector <float> > const& pente, int x){
	//cout << "x: "<<x<<endl;	
	if (x<=bpt[t][0]) return 0.0;	
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

float g_ki(int i,int k,vector< vector<float> > const& f_ks,vector<vector <float> > const& pente, int tau,int bk){
	return pente[k][i]*tau+func_fks(f_ks,k+1,tau,bk);
	
}

/*
float lotsizdyn(data const& d, int choix){//vector< vector<float> > f_ki,vector< vector<float> > g_ki){//vector<int> q0){
	

	//dtToInt(d);
	vector<int> dt = dtToInt(d,choix);

	auto start = chrono::steady_clock::now();

	vector< vector<float> > f_ki (d.cardT, vector<float> (d.nb_bp[0],0.0));
	for(int k=0; k<d.cardT; ++k){
		for(int i=0; i<d.nb_bp[k]; ++i){
			f_ki[k][i] = d.valbpt[k][i] - d.pente[k][i]*d.bpt[k][i];
			//cout << f_ki[k][i] <<endl;
		}
	}

	//vector<int> bk (d.cardT, batMax);
	vector< vector<float> > F_ks (d.cardT, vector<float> (batMax+1,0.0));
	vector< vector<int> > x_ks (d.cardT, vector<int> (batMax+1,0));

	//calcul F_Ts
	//cout <<"F["<<d.cardT-1<<"] = [";
	for (int s=0; s<=batMax; ++s){
		if(s<=dt[d.cardT-1]) F_ks[d.cardT-1][s] = p_t(d.cardT-1,d.nb_bp,d.bpt,d.valbpt,d.pente,dt[d.cardT-1]-s);
		else F_ks[d.cardT-1][s] = infini_float;
		//cout << F_ks[d.cardT-1][s] <<" ";
	}
	//cout << "]"<<endl;
	int lastk=d.cardT-1;
	int k = d.cardT-1;
	//deque<vector< vector<float> > > G_ki;
	//deque<vector<vector<int> > > x_kis;
	while ((k!=0)&&(chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count()/1000000.0) < 300.0){
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
			while((s<=batMax)&&(seq!=true)){
				//cout <<"while s : " << s << endl;
				int lb = max((int)d.bpt[k][i]+1+s-dt[k],0);
				//int lb = max(d.bpt[k][i]+s-dt[k],0);
				//int ub = min(d.bpt[k][i+1]-1+s-dt[k],batMax);
				int ub = min((int)d.bpt[k][i+1]+s-dt[k],batMax);
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
						if(g_ki(i,k,F_ks,d.pente,ub,batMax)<g_ki(i,k,F_ks,d.pente,u[0],batMax)){
							u.push_front(ub);
							--ub;
						}
						else{
							--ub;
						}
					} // ici, les u sont triés selon g(u1) < g(u2) < ...
											
					Gki.push_back(g_ki(i,k,F_ks,d.pente,u[0],batMax));
					if(dt[k]-s+u[0]>=0) xki_s.push_back(dt[k]-s+u[0]);
					else xki_s.push_back(infini_int);
					++s;
					while(s<=batMax){ // on parcours les s et on met à jour la sequence u
							
						//if((u.size()!=0)&&((d.bpt[k][i]+1+s-dt[k])==u[0])){
						if((u.size()!=0)&&((d.bpt[k][i]+s-dt[k])==u[0])){
							u.pop_front();
						}
						if((d.bpt[k][i+1]+s-dt[k])<=batMax){
							while((u.size()!=0)&&(g_ki(i,k,F_ks,d.pente,u.back(),batMax)>=g_ki(i,k,F_ks,d.pente,d.bpt[k][i+1]+s-dt[k],batMax))){
								u.pop_back();
							}
							u.push_back(d.bpt[k][i+1]+s-dt[k]);
						}
						//cout << "SEQUENCE U MAJ NUMERO "<<s<<endl;
						//for(int j=0;j<u.size();++j){ cout << u[j] << " ";}
						//cout << endl;
						if(u.size()!=0){
							Gki.push_back(g_ki(i,k,F_ks,d.pente,u[0],batMax));
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
		for(int s=0;s<=batMax;++s){
			float fGauche;
			fGauche = func_fks(F_ks,k+1,s-dt[k],batMax);
			float minFdroite = f_ki[k][0]+d.pente[k][0]*(dt[k]-s)+Gk_i[0][s];
			int imin = 0;			
			for(int i=1;i<d.nb_bp[k];++i){
				if((f_ki[k][i]+d.pente[k][i]*(dt[k]-s)+Gk_i[i][s]>=0)&&(f_ki[k][i]+d.pente[k][i]*(dt[k]-s)+Gk_i[i][s]<minFdroite)){
					minFdroite = f_ki[k][i]+d.pente[k][i]*(dt[k]-s)+Gk_i[i][s];
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
		lastk = k;
		
	}
	if((chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count()/1000000.0) < 300.0){
		/*for(int k=0; k<d.cardT; ++k){
			for(int s=0; s<=batMax; ++s){
				cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
			}
		}
		
		//cout << "F[0][0] = " <<F_ks[0][0] <<endl;
		//'''determination de xt'''    CALCUL DES XT ET DU STOCK BATTERIE
		vector<int> rt (d.cardT,0);
		vector<int> stock (d.cardT,0);// = q0;
		vector<int> xt (d.cardT,0);
		xt[0] = x_ks[0][0];
		//cout << "xt0 : " << xt[0] << endl;
		stock[0] = xt[0] - dt[0];
		//cout << "stock" << stock[0] << endl;
		xt[1] = x_ks[1][stock[0]];
		//cout << "test" << endl;
		for (int z=1; z<d.cardT-2; ++z){	
			stock[z]=(stock[z-1]+xt[z]-dt[z]);
			xt[z+1] = (x_ks[z+1][stock[z]]);
			//cout << "stock" << z << " : " << stock[z] << endl;
			//cout << "xt" << z << " : " << xt[z+1] << endl;
			
		}
		stock[d.cardT-2] = (stock[d.cardT-3]+xt[d.cardT-2]-dt[d.cardT-2]);
		xt[d.cardT-1] = dt[d.cardT-1] - stock[d.cardT-2];
		stock[d.cardT-1] = (stock[d.cardT-2]+xt[d.cardT-1]-dt[d.cardT-1]);

		rt[0] = stock[0] - d.s0;
		cout << "rt" << endl;
		cout << rt[0] << endl;
		for(int t=1; t<d.cardT; ++t){
			rt[t] = stock[t] - stock[t-1];
			//cout << rt[t] << endl;
		}

		/*for(int i=0;i<d.cardT;++i){
			cout << "demande, production, stock" << i << " : " << dt[i] << " " << xt[i] << " "<<stock[i] << endl;
		}
		for(int i=0;i<d.cardT;++i){
			cout << "xt" << i << " : " << xt[i] << endl;	
		}*/
		/*float cout_tot = 0;
		for(int i=0;i<d.cardT;++i){
			cout_tot += cout_reel(i,d.nb_bp,d.bpt,d.vald.bpt,d.pente,xt[i]);
		}
		//cout << "cout total : "<<cout_tot<<endl;
		//cout << "obj lotsizdyn : "<<F_ks[0][0]<<endl;
		return F_ks[0][0];
	} else return (float) lastk;
}*/














float lotsizdyn(data const& d, param const& p, vector<float> & var, float &tps,const vector<float> &demande){//vector< vector<float> > f_ki,vector< vector<float> > g_ki){//vector<int> q0){
	//debut timer
	
		


	vector<int> dt = dtToInt(d,demande);
	int batMax;
    //if(p.choix_dt_ls == 2) batMax = d.Q*100;
    batMax = d.Q;
	
	auto start_time = chrono::steady_clock::now();
	/*for(int i=0; i<d.cardT; ++i){
		cout << " dt["<<i<<"] = " << dt[i];
	}
	cout << endl;*/

	//auto start = chrono::steady_clock::now();

	vector< vector<float> > f_ki (d.cardT, vector<float> (d.nb_bp[0],0.0));
	for(int k=0; k<d.cardT; ++k){
		for(int i=0; i<d.nb_bp[k]; ++i){
			f_ki[k][i] = d.valbpt[k][i] - d.pente[k][i]*d.bpt[k][i];
			//cout << f_ki[k][i] << " ";
		}
		//cout << endl;
	}

	//vector<int> bk (d.cardT, batMax);
	vector< vector<float> > F_ks (d.cardT, vector<float> (batMax+1,0.0));
	vector< vector<int> > x_ks (d.cardT, vector<int> (batMax+1,0));

	//calcul F_Ts
	//cout <<"F["<<d.cardT-1<<"] = [";
	for (int s=0; s<=batMax; ++s){
		if(s<=dt[d.cardT-1]) F_ks[d.cardT-1][s] = p_t(d.cardT-1,d.nb_bp,d.bpt,d.valbpt,d.pente,dt[d.cardT-1]-s);
		else F_ks[d.cardT-1][s] = infini_float;
		//cout << F_ks[d.cardT-1][s] <<" ";
	}
	//cout << "]"<<endl;
	//int lastk=d.cardT-1;
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
			while((s<=batMax)&&(seq!=true)){
				//cout <<"while s : " << s << endl;
				int lb = max((int)d.bpt[k][i]+1+s-dt[k],0);
				//int lb = max(d.bpt[k][i]+s-dt[k],0);
				//int ub = min(d.bpt[k][i+1]-1+s-dt[k],batMax);
				int ub = min((int)d.bpt[k][i+1]+s-dt[k],batMax);
				//cout << "lb, ub : "<<lb <<" "<<ub<<endl;
				if(lb>ub){
					//cout << "SEQUENCE IMPOSSIBLE NUMERO "<<s<<endl;
					//cout <<lb<<" "<<ub<<endl;
					Gki.push_back(infini_float);
					++s;					
					xki_s.push_back(infini_int);
					//xki_s.push_back(0);
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
						if(g_ki(i,k,F_ks,d.pente,ub,batMax)<g_ki(i,k,F_ks,d.pente,u[0],batMax)){
							u.push_front(ub);
							--ub;
						}
						else{
							--ub;
						}
					} // ici, les u sont triés selon g(u1) < g(u2) < ...
											
					Gki.push_back(g_ki(i,k,F_ks,d.pente,u[0],batMax));
					if(dt[k]-s+u[0]>=0) xki_s.push_back(dt[k]-s+u[0]);
					else xki_s.push_back(infini_int);
					//else xki_s.push_back(0);
					++s;
					while(s<=batMax){ // on parcours les s et on met à jour la sequence u
							
						//if((u.size()!=0)&&((d.bpt[k][i]+1+s-dt[k])==u[0])){
						if((u.size()!=0)&&((d.bpt[k][i]+s-dt[k])==u[0])){
							u.pop_front();
						}
						if((d.bpt[k][i+1]+s-dt[k])<=batMax){
							while((u.size()!=0)&&(g_ki(i,k,F_ks,d.pente,u.back(),batMax)>=g_ki(i,k,F_ks,d.pente,d.bpt[k][i+1]+s-dt[k],batMax))){
								u.pop_back();
							}
							u.push_back(d.bpt[k][i+1]+s-dt[k]);
						}
						//cout << "SEQUENCE U MAJ NUMERO "<<s<<endl;
						//for(int j=0;j<u.size();++j){ cout << u[j] << " ";}
						//cout << endl;
						if(u.size()!=0){
							Gki.push_back(g_ki(i,k,F_ks,d.pente,u[0],batMax));
							if(dt[k]-s+u[0]>=0) xki_s.push_back(dt[k]-s+u[0]);
							else xki_s.push_back(infini_int);
							//else xki_s.push_back(0);
						}
						else{
							Gki.push_back(infini_float);
							xki_s.push_back(infini_int);
							//xki_s.push_back(0);
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
		for(int s=0;s<=batMax;++s){
			float fGauche;
			fGauche = func_fks(F_ks,k+1,s-dt[k],batMax);
			float minFdroite = f_ki[k][0]+d.pente[k][0]*(dt[k]-s)+Gk_i[0][s];
			int imin = 0;			
			for(int i=1;i<d.nb_bp[k];++i){
				if((f_ki[k][i]+d.pente[k][i]*(dt[k]-s)+Gk_i[i][s]>=0)&&(f_ki[k][i]+d.pente[k][i]*(dt[k]-s)+Gk_i[i][s]<minFdroite)){
					minFdroite = f_ki[k][i]+d.pente[k][i]*(dt[k]-s)+Gk_i[i][s];
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
		//lastk = k;
		
	}
	//if((chrono::duration_cast<chrono::microseconds>(chrono::steady_clock::now() - start).count()/1000000.0) < 300.0){
		/*for(int k=0; k<d.cardT; ++k){
			for(int s=0; s<=batMax; ++s){
				cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
			}
		}*/
		
		//cout << "F[0][0] = " << F_ks.at(0).at(0) <<endl;
		//'''determination de xt'''    CALCUL DES XT ET DU STOCK BATTERIE
		vector<float> rt (d.cardT,0.0);
		vector<int> stock (d.cardT,0);// = q0;
		vector<int> xt (d.cardT,0);
		xt[0] = x_ks[0][0];
		//cout << "xks0 : " << x_ks[0][0] << endl;
		stock[0] = xt[0] - dt[0];
		//cout << "stock" << stock[0] << endl;
		xt[1] = x_ks[1][stock[0]];
		//cout << "test" << endl;
		for (int z=1; z<d.cardT-2; ++z){	
			stock[z]=(stock[z-1]+xt[z]-dt[z]);
			xt[z+1] = (x_ks[z+1][stock[z]]);
			//cout << "stock" << z << " : " << stock[z] << endl;
			//cout << "xt" << z << " : " << xt[z+1] << endl;
			
		}
		stock[d.cardT-2] = (stock[d.cardT-3]+xt[d.cardT-2]-dt[d.cardT-2]);
		xt[d.cardT-1] = dt[d.cardT-1] - stock[d.cardT-2];
		stock[d.cardT-1] = (stock[d.cardT-2]+xt[d.cardT-1]-dt[d.cardT-1]);


		auto end_time = chrono::steady_clock::now();
		tps = chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

		//if(p.choix_dt_ls != 2){
			rt[0] = stock[0] - d.s0;
			//cout << "rt" << endl;
			//cout << rt[0] << endl;
			for(int t=1; t<d.cardT; ++t){
				rt[t] = stock[t] - stock[t-1];
				//cout << rt[t] << endl;
			}
		/*}
		else{
			rt[0] = (stock[0] - d.s0)/100;
			//cout << "rt" << endl;
			//cout << rt[0] << endl;
			for(int t=1; t<d.cardT; ++t){
				rt[t] = (stock[t] - stock[t-1])/100;
				//cout << rt[t] << endl;
			}
		}*/

		/*for(int i=0;i<d.cardT;++i){
			cout << "demande, production, stock" << i << " : " << dt[i] << " " << xt[i] << " "<<stock[i] << endl;
		}
		for(int i=0;i<d.cardT;++i){
			cout << "xt" << i << " : " << xt[i] << endl;	
		}*/
		/*float cout_tot = 0;
		for(int i=0;i<d.cardT;++i){
			cout_tot += cout_reel(i,d.nb_bp,d.bpt,d.vald.bpt,d.pente,xt[i]);
		}
		//cout << "cout total : "<<cout_tot<<endl;
		//cout << "obj lotsizdyn : "<<F_ks[0][0]<<endl;
		return F_ks[0][0];
	//} else return (float) lastk;*/
	var = rt;
	return F_ks[0][0];

}

/*
int main(){
	//creation donnees test
	int cardT = 10; 
	int qmax = 20;
	vector<int> dt = {5,10,15,20,25,12,45,24,12,65};
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
