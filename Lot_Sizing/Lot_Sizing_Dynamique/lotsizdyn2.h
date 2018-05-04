#ifndef LOTSIZDYN_H
#define LOTSIZDYN_H

#include <vector>
#include "data_struct.h"

using namespace std;

float func_fks(vector< vector<float> > const& f_ks,int k,int s,int bk);
float p_t(int t,vector<int> const& nb_bp,vector<vector <int> > const& bpt,vector< vector<float> > const& valbpt,vector<vector <float> > const& pente, int x);
float g_ki(int i,int k,vector< vector<float> > const& f_ks,vector<vector <float> > const& pente, int tau,int bk);
//void lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector< vector<float> > valbpt,vector< vector<int> > bpt,vector<vector <float> > pente,int q_max);
float lotsizdyn(data const& d,param const& p,vector<float> &var, float &tps);
//float lotsizdyn(data const& d,int choix);


#endif
