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
	int cardT, cardJ, cardM, cardR, Q;
	int s0 = 0;
	vector<int> pj;
	//vector<vector<float> > cjr;
	//vector<vector<float> > Ckr;
	//vector<float> Dk;
	vector<float> Dj;
	//vector<vector<float> > Dj;
	vector<int> rj;
	vector<int> dj;
	float consoTot;
	int releaseDateMin;
};

struct param{
	int ajout_breakpoint_voisin = 0; // 0 par defaut, 1 si on ajoute le breakpoint voisin de la colonne qui est en train d'etre ajoutee
	int ensemble_multiple = 1; // 1 par defaut, 1 si on continue de rechercher des colonnes a ajouter apres un ajout, 0 si on s'arrete
	int taches_avec_fenetre_temps = 1; // 1 si les taches ont des fenetres de temps a respecter
	//int boucle_relaxation; // 1 si on utilise la relaxation pour premiere solution ordo
	int aff_log_ordo = 0; // 1 si on souhaite afficher le deroulement cplex ordo
	//int aff_log_lotsizingcont_cplex; // 1 si on souhaite afficher le deroulement cplex lotsizcont
	//int aff_log_lotsizingent_cplex; // 1 si on souhaite afficher le deroulement cplex lotsizent
	int aff_log_compact = 0; // 1 si on souhaite afficher le deroulement cplex du modele compact
	int nb_threads_cplex = 1; // nombre de thread cplex
	int time_limit_ordo = 60; // limite de temps ordo plex en secondes
	//int time_limit_lotsizingcont; // limite de temps lotsizingcont plex en secondes
	//int time_limit_lotsizingent; // limite de temps lotsizingent plex en secondes 
	//int time_limit_lotsizingdyn; // limite de temps lotsizingdyn plex en secondes
	int time_limit_compact = 7200;
	int time_limit_gencol = 7200;
	//int time_limit_relax;
	//int nb_iter_max_boucle;
	//int choix_dt_ls;
	int type_cap = 2; //1 si la capacite est >= 0, 2 si elle est comprise entre [0, 1] (par ex 0,2 = 20%)  
	float val_cap = 0.1; //par defaut, la capacite de la batterie vaut 10% de la demande totale en energie
	//float qmin;
	//float qmax;
	//float qinit;
};

void modifPWL(data& d, vector<int> const& rt);
void modifPWL(data& d, vector<float> const& rt);
vector<int> dtToInt(data const& d, const vector<float> &demande);
float roundd(double var, int nbdec);
void initPWD(data const& init, data &d);
void initCap(data& d, param const& p);

#endif
