#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
//#include <limits>
//#include <math.h>
//#include <algorithm>
#include <chrono>
//#include <deque>
#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include <ilcplex/ilocplex.h>
#include "../struct.h"
#include "struct_gencol.h"
#include "pricer.h"

using namespace std;

#define PRICER_NAMESP1         "generatetestnewobjects"
#define PRICER_DESC            "objects generator for colgen on scheduling pb"
#define PRICER_PRIORITY        5000000
#define PRICER_DELAY           FALSE     /* only call pricer if all problem variables have non-negative reduced costs */
// pricer_delay mis à false -> notre pricer est appele si le pricer de scip n'a pas trouve de colonnes
bool verifSol(structGenCol const& sGC){
    for(int j=0; j<sGC.d.cardJ; ++j){
        int debut=0;
        double duree=0.0;
        int fin;
        //while(SCIPisEQ(sGC.scip,SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][debut]),0)){ ++debut; }
        while(sGC.varX_it[j][debut] == NULL){ ++debut; }
        if(debut<sGC.d.rj[j]){ cout<<"tache "<<j<<" debut "<<debut<<"<"<<sGC.d.rj[j]<<endl; return false; }
        for(int k=debut; k<sGC.d.cardT; ++k){
            if(sGC.varX_it[j][k] != NULL){
                if(SCIPisEQ(sGC.scip,SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][k]),1)) fin = k;
                duree += SCIPgetSolVal(sGC.scip,sGC.sol,sGC.varX_it[j][k]);
            }
            //cout << "duree "<<duree<<endl;
        }
        if(fin>=sGC.d.dj[j]){ cout<<"fin depassée pour tache "<<j<<endl; return false; }
        if(SCIPisLT(sGC.scip,duree,sGC.d.pj[j])){ cout<<"tache "<<j<<"  duree "<<duree<<"<"<<sGC.d.pj[j]<<endl; return false; }
    }
    return true;
}

double verifCoutReduit(structGenCol const& sGC, IloNumArray ui, int t){
    double cr = 0.0;
    float sumi = 0.0;
    for(int i=0; i<ui.getSize(); ++i){
        if(SCIPisEQ(sGC.scip, ui[i], 1.0)){
            cr += (sGC.u_i[i] - sGC.w_it[i][t-sGC.d.rj[i]]);
            sumi += sGC.d.Dj[i];
        }
    }
    cr += -sGC.v_t[t] - p_t(t,sumi,sGC);
    return -cr;
}

double verifCoutReduit(structGenCol const& sGC, vector<int> ui, int t){
    double cr = 0.0;
    float sumi = 0.0;
    for(auto const& i : ui){
        cr += (sGC.u_i[i] - sGC.w_it[i][t-sGC.d.rj[i]]);
        sumi += sGC.d.Dj[i];
    }
    cr += -sGC.v_t[t] - p_t(t,sumi,sGC);
    return -cr;
}

/*! \brief Method generating the new column found by the subproblem. It creates a new instance of FeasibleSet object
 This version only adds the column for the specified time
 It modifies the constraints in the restricted problem linked to the new column addition
 *
 \param[in,out] structGen Reference to structGenCol the entire problem data structure
 \param[in] valX listing alpha variable values found by solving the subproblem with Cplex
 \param[in] valY listing beta variable values found by solving the subproblem with Cplex
 \param[in] valZ listing gamma variable values found by solving the subproblem with Cplex
 \param[in] objvalue : double variable corresponding to the objective functon value found with Cplex
 *
 * \return status: a variable of type SCIP_RETCODE, equal to SCIP_OKAY if everything went well
 */
SCIP_RETCODE addObjectColumnInModel (structGenCol &sGC,IloNumArray const& valUi,int timedd,double objvalue){
	SCIP_RETCODE status = SCIP_OKAY; // devrait etre modifie lors de l ajout d une colonne, sinon erreur
	//vector<Task> listeTask;
	vector<int> taskList;
    feasibleSet l;
    bool stop=false;
    
	// Creation de la nouvelle colonne
	for(int i = 0 ; i < sGC.d.cardJ ; i++){
		//cout << "tache i=" << i << " valX="<< valUi[i] << endl;
		if( SCIPisEQ(sGC.scip, valUi[i], 1.0) ){
            //cout << "ok" << endl;
            taskList.push_back(i);
            if(taskList.size() == 1){
                l.energyDemand = sGC.d.Dj[i];
                l.deadLine = sGC.d.dj[i];
                l.releaseTime = sGC.d.rj[i];
            }
            else{
                l.energyDemand += sGC.d.Dj[i];
                if(l.deadLine > sGC.d.dj[i]) l.deadLine = sGC.d.dj[i];
                if(l.releaseTime < sGC.d.rj[i]) l.releaseTime = sGC.d.rj[i];
            }
        }
    }

    if(taskList.size() > 0){
        l.id = sGC.L.size();
        l.tasksList = taskList;
        l.cost = p_t(timedd,l.energyDemand,sGC);
        //cout << "l, cout : " << l.id <<"  "<<l.cost << endl;
        int testPres = checkSet(l,sGC);
        int id;
        bool createOk = false;
        if(testPres == -1){
            sGC.L.push_back(l);
            sGC.cardL = sGC.L.size();
            addA_il(l,sGC);
            addL_t(l,sGC);
            id = sGC.L.size()-1;
            createOk = true;
        }
        else{
            id = testPres;
            if(sGC.varY_lt[id][timedd] == NULL) createOk = true;
        }
        if(createOk){
            //vector<double> valduales;
            //vector<double> valdualesrecup;
            //vector<double> coefduales;
            // Creation de la nouvelle variable pour le moment t obtenu a l issu du sous probleme
            string name = "y_lt"+to_string(id)+","+to_string(timedd);
            //cout <<"ajout variable "+name<<endl;
            SCIP_VAR * var;
            SCIPcreateVarBasic(sGC.scip, &var, name.c_str() ,0,SCIPinfinity(sGC.scip),sGC.L[id].cost,SCIP_VARTYPE_CONTINUOUS);
            for(const auto& j : taskList){
                //Mise a jour cons1
                //cout<<"j : "<<j<<" timedd : "<<sGC.timedd<<" rj : "<<sGC.d.rj[j] << " dj : "endl;
                //SCIPprintCons(sGC.scip,sGC.cons_1[j][sGC.time-sGC.d.rj[j]],NULL);
                SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_1[j][timedd-sGC.d.rj[j]],var,-1));
                //valduales.push_back(SCIPgetDualsolLinear(sGC.scip,sGC.cons_1[j][timedd-sGC.d.rj[j]]));
                //valdualesrecup.push_back(sGC.w_it[j][timedd-sGC.d.rj[j]]);
                //coefduales.push_back(-1);
                // Mise a jour cons2
                SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_2[j],var,1));
                //valduales.push_back(SCIPgetDualsolLinear(sGC.scip,sGC.cons_2[j]));
                //valdualesrecup.push_back(sGC.u_i[j]);
                //coefduales.push_back(1);
            }
            // Mise a jour cons3
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_3[timedd],var,-1));
            //valduales.push_back(SCIPgetDualsolLinear(sGC.scip,sGC.cons_3[timedd]));
            //valdualesrecup.push_back(sGC.v_t[timedd]);
            //coefduales.push_back(-1);
            SCIP_CALL(status = SCIPaddPricedVar(sGC.scip, var, 1.0));
            /*cout << "cout reduit scip : " << SCIPgetVarRedcost(sGC.scip,var) << endl;
            cout << "cout reduit 'a la main' : " << verifCoutReduit(sGC,taskList,timedd) << endl;
            if(!SCIPisEQ(sGC.scip,objvalue,SCIPgetVarRedcost(sGC.scip,var))){
                cout << "cost : " << sGC.L[id].cost << endl;
                for(int v=0; v<valduales.size(); ++v) cout << v << " :     valduale : " << valduales[v] <<" valdualerecup : " << valdualesrecup[v] << "   coeff  : " << coefduales[v] << endl;
            }*/
            //assert(SCIPisEQ(sGC.scip,objvalue,SCIPgetVarRedcost(sGC.scip,var)));
            if(testPres==-1){
                vector<SCIP_VAR*> tab (sGC.d.cardT, NULL);
                tab[timedd] = var;
                sGC.varY_lt.push_back(tab);
            }
            else{
                sGC.varY_lt[id][timedd] = var;
            }
            sGC.nbcolgenerated++;
        }
    }
    return status;
}


/*! \brief Method creating and solving the subproblem for a predefined t
 It gives the subset of tasks in execution at time t (use a piece-wise linear obj function for the SP)
Then the column (subset of tasks) is added to the restricted master problem with one of the three adding methods
 *
 \param[in,out] pbdata Reference to structGenCol the entire problem data structure
 *
 * \return status: a variable of type SCIP_RESULT, equal to SCIP_SUCCESS if everything went well
 */
SCIP_RESULT Pr_SP1(structGenCol &sGC){
	
    int timedd = 0;
    double objfctvalue = 0.0;
    int statusint = 0;
    bool colgen = false;
    bool stop = false;
    SCIP_RESULT status = SCIP_SUCCESS;
	SCIP_RETCODE status1 = SCIP_ERROR;
    int decal; // 1 -> d.bpt[i][0] > 0   , 2 -> sinon
    
    while((timedd<sGC.d.cardT)&&(!stop)){

        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);
        if(sGC.p.aff_log_ordo == 0) cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::Threads,1);

        //donnees pwl

        /*for(int i=0; i<d.cardT; ++i){
            if(d.bpt[i][0] > 0){
                IloNumArray pente(env, d.nb_bp[i]+1);
                IloNumArray bpt(env, d.nb_bp[i]);
                pente[0] = 0.0;
                bpt[0] = d.bpt[i][0];
                for(int j=1; j<d.nb_bp[i]+1;++j){
                    pente[j] = d.pente[i][j-1];
                    if (j<d.nb_bp[i]) bpt[j] = d.bpt[i][j];
                }
                decal = 1;
                //model.add(ct[i] == IloPiecewiseLinear(xt[i],bpt,pente,0, 0));
            }
            else{
                IloNumArray pente(env, d.nb_bp[i]);
                IloNumArray bpt(env, d.nb_bp[i]-1);
                for(int j=0; j<d.nb_bp[i];++j){
                    pente[j] = d.pente[i][j];
                    if(j<d.nb_bp[i]-1) bpt[j] = d.bpt[i][j+1];
                }
                decal = 2;
                //model.add(ct[i] == IloPiecewiseLinear(xt[i],bpt,pente,d.bpt[i][0], d.valbpt[i][0]));
            }
	    }*/
        IloNumArray bpt(env, sGC.d.nb_bp[timedd]-1);
        for(int i=0; i<(sGC.d.nb_bp[timedd]-1);++i){
            bpt[i] = sGC.d.bpt[timedd][i+1];
        }

        IloNumArray pente(env, sGC.d.nb_bp[timedd]);
        for(int i=0; i<sGC.d.nb_bp[timedd];++i){
            pente[i] = sGC.d.pente[timedd][i];
        }

        //ajout variables u_i
        IloNumVarArray alpha_i (env,sGC.d.cardJ,0.0,1.0,ILOBOOL);
        
        IloNumArray valAlphai(env);

        //contrainte sum_i bi*ui - sum_k o_k2(pk+1)*vk <= qmax
        IloExpr sum_i1(env);
        IloExpr sum_i2(env);
        for(int i=0; i<sGC.d.cardJ;++i){
            if((sGC.d.rj[i]<=timedd)&&(timedd<sGC.d.dj[i])){
                sum_i1 += alpha_i[i] * (sGC.u_i[i] - sGC.w_it[i][timedd-sGC.d.rj[i]]);
                sum_i2 += alpha_i[i] * sGC.d.Dj[i];
            }
            else model.add(alpha_i[i] <= 0);
        }



        IloExpr obj(env);
        obj = sum_i1 - sGC.v_t[timedd] - IloPiecewiseLinear(sum_i2,bpt,pente,sGC.d.bpt[timedd][0],0.0);
        //if(decal==1) obj = sum_i1 - sGC.v_t[timedd] - IloPiecewiseLinear(sum_i2,bpt,pente,0.0,0.0);
        //else obj = sum_i1 - sGC.v_t[timedd] - IloPiecewiseLinear(sum_i2,bpt,pented.bpt[i][0], d.valbpt[i][0]);
        model.add(IloMaximize(env, obj));
        obj.end();
        
			
        // cout << "after model creation to delete" << endl;

        auto start_time = chrono::steady_clock::now();
        if (!cplex.solve()){
            cout << "Failed to optimize LP." << endl;
            statusint = -1;
            cout << cplex.getStatus() << endl;
            cout << cplex.getCplexStatus() << endl;
            cout <<	cplex.getCplexSubStatus() << endl;
        }

        // cout << "after model solve to delete" << endl;

        auto end_time = chrono::steady_clock::now();
        sGC.tpsSP += chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;
        if(cplex.getStatus() != IloAlgorithm::Optimal){
            cerr << "Cplex did not find the optimal solution for the subproblem !!!!!!!!!! " << endl;
            cout << "elapsed time=" << SCIPgetSolvingTime(sGC.scip) << endl;
            cout << cplex.getStatus() << endl;
            cout << cplex.getCplexStatus() << endl;
            cout <<	cplex.getCplexSubStatus() << endl;
            cout << "obj =" << cplex.getObjValue() << endl;
            statusint = -1;
        }
        else{
            objfctvalue = cplex.getObjValue();
            cplex.getValues(valAlphai,alpha_i);
            //if(pbdata.Params.CHOIX_AFFICHAGE==1){
                //cout << "Solution status = " << cplex.getStatus() << endl;
                //cout << "Solution value = " << cplex.getObjValue() << endl;
                //cout << "Values alpha = " << valX << endl;
            //}
        }
        if(statusint==0){
            if (SCIPisGT(sGC.scip,objfctvalue,0))
            {
                //cout << "temps export : " << timedd << endl;
                //cplex.exportModel("SPNoStock.lp");
                //inittime = sGC.time;
                // set result pointer
                colgen=true;
                /*for(int i=0; i<valUi.getSize(); i++){
                    cout << "u["<<i<<"] = "<<valUi[i] << endl;
                }
                for(int k=0; k<valVk.getSize(); k++){
                    cout << "v["<<k<<"] = "<<valVk[k] << endl;
                }
                cout<<"time = "<<sGC.time<<endl;*/
                status = SCIP_SUCCESS;
                //cout << "cout reduit cplex : " << -objfctvalue << endl;
                //cout << "cout reduit SP CPLEX : " << objfctvalue << endl;
                status1 = addObjectColumnInModel(sGC,valAlphai,timedd,-objfctvalue);
                assert(status1==1);
                
                timedd++;
                //assert(pbdata.time - inittime >= 0);
                //cpteur = cpteur + (pbdata.time - inittime);
                if(sGC.p.ensemble_multiple==0) stop = true;
            }else
            {   
                //cout << "Cout reduit non négatif au temps "<<timedd<< endl;
                status = SCIP_SUCCESS;
                timedd++;
                //cpteur = cpteur + 1;
            }
        }
        else{
            if(!colgen) status = SCIP_DIDNOTRUN;
            else status = SCIP_SUCCESS;
            stop = true;
        }
        env.end();
    }
    
	return status;
}





/*! \typedef SCIP_DECL_PRICERREDCOST
 \brief Reduced cost pricing method of variable pricer for feasible LPs.
 *
 */
static
SCIP_DECL_PRICERREDCOST(pricerRedcost)
{
	structGenCol* pbdata; /* the data of the pricer */
	SCIP_PROBDATA* probdata;
	SCIP_RETCODE status = SCIP_OKAY;
	assert(scip != NULL);
	assert(pricer != NULL);
	probdata = SCIPgetProbData(scip);
	pbdata = (structGenCol*) probdata;
	auto start_time = chrono::steady_clock::now();
	//int duration = pbdata->instance.getDeadLine()-pbdata->instance.getReleaseTime();

	/*if(pbdata->Params.ACTIVATE_VERIF==1){
		assert(Pr_VerifSol(*pbdata,duration)==true);
	}*/
	pbdata->cptPricer++;
    //cout << " nb pricer vars : "<< SCIPgetNPricevarsFound(scip) << "   nbcolgen :" << pbdata->nbcolgenerated << endl;
	// --------------------------------------   recuperation des variables duales
    int fois = 0;
    for(int i=0; i<pbdata->d.cardJ; ++i){
        pbdata->u_i[i] = SCIPgetDualsolLinear(scip,pbdata->cons_2[i]);
		for(int t=0; t<pbdata->d.cardT; ++t){
            if((pbdata->d.rj[i]<=t)&&(t<pbdata->d.dj[i])){
                pbdata->w_it[i][t-pbdata->d.rj[i]] = SCIPgetDualsolLinear(scip,pbdata->cons_1[i][t-pbdata->d.rj[i]]);
            }
            if(fois==0){
                pbdata->v_t[t] = SCIPgetDualsolLinear(scip,pbdata->cons_3[t]);
            }
        }
        ++fois;
	}
	/*
	int aff = 1;
	if(aff==1){
		cout<<"lowerbound = "<<*lowerbound<<endl;
		for(int i = 0 ; i < pbdata->d.cardJ;i++)
			cout << pbdata->beta_i[i] <<" ";
		cout << endl;
		for(int i = 0 ; i < duration;i++)
			cout << pbdata->v[i] <<" ";
		cout << endl;
		for(int i = 0 ; i < pbdata->instance.getNbOfTasks(); i++){
			for(int j = pbdata->instance.getTasksList().at(i).getReleaseTime() ; j < pbdata->instance.getTasksList().at(i).getDeadLine() ;j++){
				cout <<"Tache "<<i+1 << " Instant "<< j<<" : "<< pbdata->w[i][j] << endl;
			}
		}
		cout << endl;
	}*/

	// --------------------------------Construction du Solver du sous probleme (pour tout t)

	// Appel de la fonction de la fonction lineaire par morceaux
    *result = Pr_SP1(*pbdata);
		
	if (*result==SCIP_DIDNOTRUN) {
                cout << "DID NOT RUN " << endl;
                SCIPinterruptSolve(pbdata->scip);	
        }
	/*if(pbdata->Params.CHOIX_AFFICHAGE==1){
		char path1[200];
		sprintf(path1, "%s%s",pbdata->Params.PATH_PROBLEM,"NodeOfColumnGen.lp");
		SCIP_CALL(status = SCIPwriteTransProblem(scip, path1, "lp", FALSE));
	}
	assert(status==1);*/
	auto end_time = chrono::steady_clock::now();
    pbdata->tpsPricer += chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	return SCIP_OKAY;
}


/*! \fn SCIP_DECL_PRICERINITSOL(pricerInitsolSP1)
    \brief solving process initialization method of variable pricer (called when branch and bound process is about to begin) here the pointers to the transformed variables and constraints are collected to ensure that the later procedures will use the right links
 */
static
SCIP_DECL_PRICERINITSOL(pricerInitsolSP)
{
	structGenCol* pbdata; // the data of the pricer
	SCIP_PROBDATA* probdata;

	assert(scip != NULL);
	assert(pricer != NULL);

	probdata = SCIPgetProbData(scip);
	pbdata = (structGenCol *) probdata;
	// get transformed variables
    for(int l=0; l<pbdata->L.size(); ++l){
        for(int t=pbdata->L[l].releaseTime; t<pbdata->L[l].deadLine; ++t){
            SCIPgetTransformedVar(scip, pbdata->varY_lt[l][t],&(pbdata->varY_lt[l][t]));
            //assert(pbdata->varY_lkt[l][k][t] != NULL);
        }
	}

	for(int i=0; i<pbdata->d.cardJ; i++){
		for(int t = pbdata->d.rj[i] ; t < pbdata->d.dj[i] ; t++){
			SCIPgetTransformedVar(scip, pbdata->varX_it[i][t],&(pbdata->varX_it[i][t]));
			//assert(pbdata->tabVarXit.at(i).at(j-pbdata->instance.getTasksList().at(i).getReleaseTime()) != NULL);
		}
	}

	// get transformed constraints
	//int duration = pbdata->instance.getDeadLine()-pbdata->instance.getReleaseTime();
	int fois = 0;

	for(int i=0; i<pbdata->d.cardJ; ++i){
        SCIPgetTransformedCons(scip,pbdata->cons_2[i],&(pbdata->cons_2[i]));
        
		for(int t=0; t<pbdata->d.cardT; ++t){
            if((pbdata->d.rj[i]<=t)&&(t<pbdata->d.dj[i])){
                SCIPgetTransformedCons(scip,pbdata->cons_1[i][t-pbdata->d.rj[i]],&(pbdata->cons_1[i][t-pbdata->d.rj[i]]));
                //assert(pbdata->tabConsNbSets[j] != NULL);
            }
            if(fois==0){
                SCIPgetTransformedCons(scip,pbdata->cons_3[t],&(pbdata->cons_3[t]));
            }
        }
        ++fois;
	}
	return SCIP_OKAY;
}




/*! \brief Method creating the pricer, activates it and includes it in SCIP
 *
 \param[in,out] structGen Reference to structGenCol the entire problem data structure
 *
 * \return status: a variable of type SCIP_RETCODE, equal to SCIP_OKAY if everything went well
 */
SCIP_RETCODE includePricer(structGenCol &sGC){
	SCIP_PRICER* pricer;
	SCIP_PRICERDATA* pricerdata;
	pricerdata = (SCIP_PRICERDATA*) &sGC;
	SCIP_RETCODE status;
	pricer = NULL;
	assert (pricerdata != NULL);
	SCIP_CALL(status =  SCIPincludePricerBasic(sGC.scip, &pricer, PRICER_NAMESP1, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY, pricerRedcost, 0, pricerdata));
	assert(pricer != NULL);
	SCIP_CALL(status = SCIPsetPricerInitsol(sGC.scip, pricer, pricerInitsolSP));
    return status;
}