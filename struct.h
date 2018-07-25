#ifndef struct_H
#define struct_H

#include <vector>

using namespace std;


struct data{
	// PWL
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
	float consoTot;
	int releaseDateMin;
};

struct param{
	int ajout_breakpoint_voisin;
	int ensemble_multiple;
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
	int time_limit_relax;
	int nb_iter_max_boucle;
	int choix_dt_ls;
	int type_cap;
	float val_cap;
	float qmin;
	float qmax;
	float qinit;
};

void modifPWL(data& d, vector<int> const& rt);
void modifPWL(data& d, vector<float> const& rt);
vector<int> dtToInt(data const& d, int choix);
float roundd(double var, int nbdec);
void initPWD(data const& init, data &d);
void initCap(data& d, param const& p);

#endif
