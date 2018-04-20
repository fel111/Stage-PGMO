#ifndef LOTSIZDYN_H
#define LOTSIZDYN_H

#include <vector>
#include "data_struct.h"

using namespace std;

float func_fks(vector< vector<float> > f_ks,int k,int s,int bk);
float p_t(int t,vector<int> nb_bp,vector<vector <int> > bpt,vector< vector<float> > valbpt,vector<vector <float> > pente, int x);
float g_ki(int i,int k,vector< vector<float> > f_ks,vector<vector <float> > pente, int tau,int bk);
//void lotsizing(int cardT,vector<int> dt,vector<int> nb_bp,vector< vector<float> > valbpt,vector< vector<int> > bpt,vector<vector <float> > pente,int q_max);
float lotsizdyn(data d,int choix,vector<float> &var);
float lotsizdyn(data d,int choix);


#endif
