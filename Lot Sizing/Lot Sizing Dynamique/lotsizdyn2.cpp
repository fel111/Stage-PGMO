#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <limits>
#include <deque>

using namespace std;

float infini_float = numeric_limits<float>::max();
int infini_int = numeric_limits<int>::max();

float func_fks(vector< vector<float> > f_ks,int k,int s,int bk){ //renvoie la valeur de f_ks
	if ((s<0) || (s>bk)) return infini_float;
	else return f_ks[k][s];
}

float p_t(int t,vector<int> nb_bp,vector<vector <int> > bpt,vector< vector<float> > valbpt,vector<vector <float> > pente, float x){
	//cout << "x: "<<x<<endl;	
	if (x==0.0) return 0.0;	
	for (int i=0; i<nb_bp[t]; ++i){
		//cout<<bpt[t][i]<<"<"<<x<<"<="<<bpt[t][i+1]<<endl;
		if((x>bpt[t][i])&&(x<=bpt[t][i+1])){
			//cout << "valreel" <<  pente[t][i]*x+valbpt[t][i]<<endl;
			return pente[t][i]*x+valbpt[t][i];
		}	
	}
}

float g_ki(int i,int k,vector< vector<float> > f_ks,vector<vector <float> > pente, int tau,int bk){
	return pente[k][i]*tau+func_fks(f_ks,k+1,tau,bk);
	
}

/*float G_ik(int i, int k, int eta, vector< vector<int> > bpt,int & u0, bool upd,vector<int> dt,int s,vector< vector<float> > valbpt,int q_max){
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

vector<int> lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector< vector<float> > valbpt,vector< vector<int> > bpt,vector<vector <float> > pente,int q_max){//vector< vector<float> > f_ki,vector< vector<float> > g_ki){//vector<int> q0){
	//debut timer
	//auto start = chrono::high_resolution_clock::now();

	vector<int> bk (cardT, q_max);
	vector< vector<float> > F_ks (cardT, vector<float> (q_max,0.0));
	vector< vector<float> > x_ks (cardT, vector<float> (q_max,0.0));

	//calcul F_Ts
	for (int s=0; s<bk[cardT-1]; ++s){
		if(s<=dt[cardT-1]) F_ks[cardT-1][s] = p_t(cardT-1,nb_bp,bpt,valbpt,pente,dt[cardT-1]-s);
		else F_ks[cardT-1][s] = infini_float;
		cout <<"F["<<cardT-1<<"]["<<s<<"] = "<<F_ks[cardT-1][s] <<endl;
	}
	
	int k = cardT;
	deque<vector< vector<float> > > G_ki;
	deque<vector<vector<float> > > x_kis;
	while (k>=0){
		--k;
		vector< vector<float> >Gk_i;
		vector< vector<float> >xk_is;
		for(int i=0; i<nb_bp[k]; ++i){
			vector<float> Gki ;
			vector<float> xki_s ;
			int s=0;
			while(s<=qmax){
				float lb = max(bpt[k][i]+1+s-dt[k],0);
				float ub = min(bpt[k][i+1]+s-dt[k],q_max);
				if(lb>ub){
					Gki.push_back(infini_float);
					xkis.push_back(-infini_float);
					++s;					
					//xks[i] = infini;
				}
				else{
					deque<int> u = {ub};
					--ub;
					while(lb<=ub){
						if(g_ik(i,k,F_ks,pente,ub,qmax)<g_ik(i,k,F_ks,pente,u[0],q_max)){
							u.push_front(ub);
							--ub;
						}
						else{
							--ub;
						}
					}
					Gki.push_back(g_ik(i,k,F_ks,pente,u[0],q_max));
					xkis.push_back(dt[k]-s+u[0]);
					++s;
					while(s<=q_max){
						if((u.size()!=0)&&((bpt[k][i]+s-dt[k])==u[0])){
							u.pop_front();	
						}
						if((bpt[k][i+1]+s-dt[k])<=q_max){
							while((u.size()!=0)&&(g_ik(i,k,F_ks,pente,u.back(),q_max)>=g_ik(i,k,F_ks,pente,bpt[k][i+1]+s-dt[k],q_max))){
								u.pop_back();
							}
							u.push_back(bpt[k][i+1]+s-dt[k]);
						}
						if(u.size()!=0){
							Gki.push_back(g_ik(i,k,F_ks,pente,u[0],q_max));
							xki_s.push_back(dt[k]-s+u[0]);
						}
						else{
							Gki.push_back(infini_float);
							xki_s.push_back(-infini_float);
						}
						++s;
					}
				}
			}
			
			Gk_i.push_back(Gki);
			xk_is.push_back(xki_s);
		}
		G_ki.push_front(Gk_i);
		x_kis.push_front(xk_is);
		for(int s=0;s<=q_max;++s){
			deque<float> list2;
			list2.push_back(func_fks(F_ks,k+1,s-dt[k],q_max));
			float min1 = list2[0];			
			for(int i=0;i<nb_bp[k];++i){
				if(valbpt[k][i]+pente[k][i]*(dt[k]-s)+Gk_i[i][s]>=0)
					list2.push_back(valbpt[k][i]+pente[k][i]*(dt[k]-s)+Gk_i[i][s]);
				else list2.push_back(infini_float);
				if(list2[i+1]<min1) min1=list2[i+1];				
			}
			F_ks[k].push_back(min1);
			if(min1==func_fks(F_ks,k+1,s-dt[k],q_max)) x_ks[k].push_back(0.0);
			else{
				deque<float> listi = list2;
				listi.pop_front();
				auto it = find(listi.begin(), listi.end(), old_name_);
				if (it == Names.end())
{
  // name not in vector
} else
{
  auto index = std::distance(Names.begin(), it);
}
			}

		}
		//	
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
			vector<float> Gki (nb_bp[k],infini_float);
			//vector<float> xks (nb_bp[k],infini);
			int u_min = 0;
			float valmin = infini_float;
			
				float lb = max(bpt[k][i]+1+s-dt[k],0);
				float ub = min(bpt[k][i+1]+s-dt[k],q_max);
				if(lb>ub){
					Gki[i] = infini_float;
					//xks[i] = infini;
				}
				else{
					Gki[i] = g_ki(i,k,F_ks,pente,lb,bk[k]);
					if(Gki[i] < valmin){
						u_min = i;
					}
				}
				//cout <<"G["<<k<<"]["<<i<<"] = "<<Gki[i] <<endl;
			}
			//cout << "k : "<<k<<"   s : "<<s<<"    u_min : "<<u_min<<endl;
			x_ks[k][s] = dt[k]-s+u_min;
			F_ks[k][s] = min(func_fks(F_ks,k+1,s-dt[k],bk[k]),valbpt[k][u_min]+pente[k][u_min]*(dt[k]-s)+Gki[u_min]);
			cout <<"F["<<k<<"]["<<s<<"] = "<<F_ks[k][s] <<endl;
			cout <<"x["<<k<<"]["<<s<<"] = "<<x_ks[k][s] <<endl;
		}
		--k;
	}
	//'''determination de xt'''    CALCUL DES XT ET DU STOCK BATTERIE
	vector<int> stock (cardT,0);// = q0;
	vector<int> xt (cardT,0);
	xt[0] = x_ks[0][0];
	stock[0] = xt[0] - dt[0];
	xt[1] = x_ks[1][stock[0]];
	for (int z=1; z<cardT-1; ++z){		
		stock[z]=(stock[z-1]+xt[z]-dt[z]);
		//cout << "stock" << z << " : " << stock[z] << endl;
		xt[z+1] = (x_ks[z+1][stock[z]]);
	}
	stock[cardT-1] = (stock[cardT-2]+xt[cardT-1]-dt[cardT-1]);
	for(int i=0;i<cardT;++i){
		cout << "stock" << i << " : " << stock[i] << endl;
	}	

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
	vector<vector<float> > valbpt;
	vector<vector<float> > pente;
	vector<int> nb_bp = {3,3,3,3,3};
	for(int i=0;i<cardT;++i){
		vector<float> temp_pente = {1.0,2.0,1.0};
		vector<int> temp_bpt = {0,10,20,infini_int};
		vector<float> temp_valbpt = {0.0,10.0,30.0};
		pente.push_back(temp_pente);
		bpt.push_back(temp_bpt);
		valbpt.push_back(temp_valbpt);
	}
	
	vector<int> xt = lotsizing(cardT,dt,nb_bp,valbpt,bpt,pente,qmax);
	for(int i=0; i<xt.size(); ++i){
		cout << "xt"<<i<<" : "<<xt[i] << endl;
	}
	return 0;	
}
