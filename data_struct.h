#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H

#include <vector>

using namespace std;


struct data{
	// PWD
	vector<int> nb_bp;
	vector<vector <int> > bpt;
	vector< vector<float> > valbpt;
	vector<vector<float> > pente;
	vector<float> dt;
	// Data Ordo
	int cardT, cardJ, cardM, cardR, s0, Q;
	vector<int> pj;
	vector<vector<float> > cjr;
	vector<vector<float> > Ckr;
	vector<float> Dk;
	vector<vector<float> > Djk;
	vector<int> rj;
	vector<int> dj;
};

void modifPWL(data& d, vector<int> rt);
void modifPWL(data& d, vector<float> rt);
vector<int> dtToInt(data d, int choix);
float roundd(double var, int nbdec);
void initPWD(data init, data &d);


#endif
