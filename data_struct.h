#ifndef DATA_STRUCT_H
#define DATA_STRUCT_H

#include <vector>

using namespace std;


struct data{
	// PWD
	vector<int> nb_bp;
	vector<vector <float> > bpt;
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

struct param{
	int taches_avec_fenetre_temps; // 1 si les taches ont des fenetres de temps a respecter
	int boucle_relaxation; // 1 si on utilise la relaxation pour premiere solution ordo
	int aff_log_ordo_cplex; // 1 si on souhaite afficher le deroulement cplex ordo
	int aff_log_lotsizingcont_cplex; // 1 si on souhaite afficher le deroulement cplex lotsizcont
	int aff_log_lotsizingent_cplex; // 1 si on souhaite afficher le deroulement cplex lotsizent
	int aff_log_compact_cplex; // 1 si on souhaite afficher le deroulement cplex du modele compact
	int nb_threads_cplex; // nombre de thread cplex
	int time_limit_ordo; // limite de temps ordo plex en secondes
	int time_limit_lotsizingcont; // limite de temps lotsizingcont plex en secondes
	int time_limit_lotsizingent; // limite de temps lotsizingent plex en secondes 
	int time_limit_lotsizingdyn; // limite de temps lotsizingdyn plex en secondes
	int time_limit_compact;
};

void modifPWL(data& d, vector<int> rt);
void modifPWL(data& d, vector<float> rt);
vector<int> dtToInt(data d, int choix);
float roundd(double var, int nbdec);
void initPWD(data init, data &d);


#endif
