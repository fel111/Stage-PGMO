#include <iostream>
//#include <fstream>
//#include <sstream>
#include <vector>
#include <string>
//#include <limits>
//#include <math.h>
//#include <algorithm>
//#include <chrono>
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
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

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
SCIP_RETCODE addObjectColumnInModel (structGenCol &sGC,IloNumArray valUi,IloNumArray valVk,int time){
	SCIP_RETCODE status = SCIP_ERROR; // devrait etre modifie lors de l ajout d une colonne, sinon erreur
	//vector<Task> listeTask;
	vector<int> taskList;
    feasibleSet l;
    bool stop=false;
    
	// Creation de la nouvelle colonne
	for(int i = 0 ; i < sGC.d.cardJ ; i++){
		//cout << "tache i=" << i << " valX="<< valX[i] << endl;
		if( SCIPisEQ(sGC.scip, valUi[i], 1.0) ){
            taskList.push_back(i);
            if(taskList.size() == 1){
                l.energyDemand = sGC.d.Djk[i][0];
                l.deadLine = sGC.d.dj[i];
                l.releaseTime = sGC.d.rj[i];
            }
            else{
                l.energyDemand += sGC.d.Djk[i][0];
                if(l.deadLine > sGC.d.dj[i]) l.deadLine = sGC.d.dj[i];
                if(l.releaseTime < sGC.d.rj[i]) l.releaseTime = sGC.d.rj[i];
            }
        }
    }
    
    int k = 0 ;
    while( SCIPisEQ(sGC.scip, valVk[k], 0.0) ) ++k;


    if(taskList.size() > 0){
        l.id = sGC.L.size();
        l.timeGen = time; 
        l.tasksList = taskList;
        int testPres = checkSet(l,sGC);
        int id;
        if(testPres == -1){
            sGC.L.push_back(l);
            sGC.cardL = sGC.L.size();
            addSetK_l(l,sGC);
            addA_il(l,sGC);
            addL_t(l,sGC);
            id = sGC.L.size()-1;
        }
        else id = testPres;
        // Creation de la nouvelle variable pour le moment t obtenu a l issu du sous probleme
        string name = "y_lkt"+to_string(id)+","+to_string(k)+","+to_string(time);
        cout <<"ajout variable "+name<<endl;
        SCIP_VAR * var;
        SCIPcreateVarBasic(sGC.scip, &var, name.c_str() ,0,SCIPinfinity(sGC.scip),sGC.d.valbpt[0][k],SCIP_VARTYPE_CONTINUOUS);
        for(const auto& j : taskList){
            //Mise a jour cons1
            //cout<<"j : "<<j<<" time : "<<sGC.time<<" rj : "<<sGC.d.rj[j] << " dj : "endl;
            //SCIPprintCons(sGC.scip,sGC.cons_1[j][sGC.time-sGC.d.rj[j]],NULL);
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_1[j][time-sGC.d.rj[j]],var,-1));
            // Mise a jour cons2
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_2[j],var,1));
        }
        // Mise a jour cons3
        SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_3[k/2][time],var,-1));
        // Mise a jour cons8
        //cout<<"contrainte 8 : l.energyDemand-sGC.d.bpt[0][k] = "<<l.energyDemand-sGC.d.bpt[0][k]<<endl;
        SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_8[time],var,l.energyDemand-sGC.d.bpt[0][k]));
        // Mise a jour cons9
        SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_9[time],var,l.energyDemand));
        SCIP_CALL(status = SCIPaddPricedVar(sGC.scip, var, 1.0));
        // ???????????? START
        cout << "cout reduit de colonne ajoutee : "<<SCIPgetVarRedcost(sGC.scip,var)<<endl;
        // ??????????? END
        if(testPres==-1){
            vector<vector<SCIP_VAR*> > tab (sGC.d.nb_bp[0], vector<SCIP_VAR*> (sGC.d.cardT, NULL));
            tab[k][time] = var;
            sGC.varY_lkt.push_back(tab);
        }
        else{
            sGC.varY_lkt[id][k][time] = var;
        }

        // PARTIE AJOUT COLONNE Y_LHT, AVEC H LE BREAKPOINT VOISIN DE K
        int k2;
        if(k%2==0) k2=k+1;
        else k2=k-1; 
        if((testPres==-1)||((testPres!=-1)&&(sGC.varY_lkt[id][k2][time]==NULL))){
            string name2 = "y_lkt"+to_string(id)+","+to_string(k2)+","+to_string(time);
            SCIP_VAR * var2;
            SCIPcreateVarBasic(sGC.scip, &var2, name2.c_str() ,0,SCIPinfinity(sGC.scip),sGC.d.valbpt[0][k2],SCIP_VARTYPE_CONTINUOUS);
            for(const auto& j : taskList){
                //Mise a jour cons1
                //cout<<"j : "<<j<<" time : "<<sGC.time<<" rj : "<<sGC.d.rj[j] << " dj : "endl;
                //SCIPprintCons(sGC.scip,sGC.cons_1[j][sGC.time-sGC.d.rj[j]],NULL);
                SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_1[j][time-sGC.d.rj[j]],var2,-1));
                // Mise a jour cons2
                SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_2[j],var2,1));
            }
            // Mise a jour cons3
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_3[k2/2][time],var2,-1));
            // Mise a jour cons8
            //cout<<"contrainte 8 : l.energyDemand-sGC.d.bpt[0][k] = "<<l.energyDemand-sGC.d.bpt[0][k]<<endl;
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_8[time],var2,l.energyDemand-sGC.d.bpt[0][k2]));
            // Mise a jour cons9
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_9[time],var2,l.energyDemand));
            SCIP_CALL(status = SCIPaddPricedVar(sGC.scip, var2, 1.0));
            // ???????????? START
            cout << "cout reduit de colonne ajoutee : "<<SCIPgetVarRedcost(sGC.scip,var2)<<endl;
        }
        return status;
    }
    /*else{ //cas pas de tâches donc variable y0_kt
        // Creation de la nouvelle variable pour le moment t obtenu a l issu du sous probleme
        string name = "y0_kt"+to_string(k)+","+to_string(time);

        SCIP_VAR * var;
        SCIPcreateVarBasic(sGC.scip, &var, name.c_str() ,0,SCIPinfinity(sGC.scip),sGC.d.valbpt[0][k],SCIP_VARTYPE_CONTINUOUS);

        // Mise a jour cons4
        //SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_4[k/2][time],var,-1));
        // Mise a jour cons8
        SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_8[time],var,-sGC.d.bpt[0][k]));
        SCIP_CALL(status = SCIPaddPricedVar(sGC.scip, var, 1.0));
        cout << "--------COLONNES AJOUTEE-------"<<endl;
        cout << "cout reduit de colonne ajoutee : "<<SCIPgetVarRedcost(sGC.scip,var)<<endl;
        sGC.varY0_kt[k][time] = var;
        //time++;
        return status;
    }*/
    return SCIP_OKAY;
        //}else cout << "l deja present " << l.id << endl;
	//}
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
	
    sGC.time = 0;
    double objfctvalue = 0.0;
    int statusint = 0;
    bool colgen = false;
    bool stop = false;
    SCIP_RESULT status = SCIP_SUCCESS;
	SCIP_RETCODE status1 = SCIP_ERROR;
    
    while((sGC.time<sGC.d.cardT)&&(!stop)){

        IloEnv env;
        IloModel model(env);
        IloCplex cplex(model);
        cplex.setOut(env.getNullStream());
        cplex.setParam(IloCplex::Threads,1);

        //ajout variables u_i
        IloNumVarArray u_i (env,sGC.d.cardJ,0.0,1.0,ILOBOOL);
        
        //ajout variables v_k
        IloNumVarArray v_k (env,sGC.d.nb_bp[0],0.0,1.0,ILOBOOL);
        
        IloNumArray valUi(env);
        IloNumArray valVk(env);

        //contrainte sum_i bi*ui - sum_k o_k2(pk+1)*vk <= qmax
        IloExpr sum_i1(env);
        IloExpr sum_i2(env);
        for(int i=0; i<sGC.d.cardJ;++i){
			//cout << i << "    " << sGC.d.rj[i] << " <= " << sGC.time<<" && "<< sGC.time<<" < " <<sGC.d.dj[i]<<endl;
            if((sGC.d.rj[i]<=sGC.time)&&(sGC.time<sGC.d.dj[i])){
                sum_i1 += sGC.d.Djk[i][0]*u_i[i]; //sum bi ui
                //???
                //cout << "alpha : "<<sGC.alpha_it[i] << endl;
                //cout << "time" << sGC.time << endl;
                sum_i2 += u_i[i] * (sGC.beta_i[i] //sum_i obj
                    - sGC.alpha_it[i][sGC.time-sGC.d.rj[i]] 
                    + sGC.d.Djk[i][0] 
                        * (sGC.delta_t[sGC.time]
                           + sGC.phi_t[sGC.time]));
            }
            else model.add(u_i[i] <= 0);
        }
        IloExpr sum_k1(env);
        /*for(int k=0; k<sGC.d.nb_bp[0]-2; ++k){
            if(k%2==0) sum_k1 += sGC.d.bpt[0][k+3]*v_k[k];
            else sum_k1 += sGC.d.bpt[0][k+2]*v_k[k];
        }*/
        IloExpr sum_k2(env);
        /*for(int k=2; k<sGC.d.nb_bp[0]; ++k){
            if(k%2==0) sum_k2 += sGC.d.bpt[0][k-1]*v_k[k];
            else sum_k2 += sGC.d.bpt[0][k-2]*v_k[k];
        }*/
        IloExpr sum_k3(env);
        IloExpr sum_k4(env);
        for(int k=0; k<sGC.d.nb_bp[0]; ++k){
			if(k%2==0){ 
                sum_k1 += sGC.d.bpt[0][k]*v_k[k];
                sum_k2 += sGC.d.bpt[0][k+1]*v_k[k];
            }
            else{
                sum_k1 += sGC.d.bpt[0][k-1]*v_k[k];
                sum_k2 += sGC.d.bpt[0][k]*v_k[k];
            }
            sum_k3 += v_k[k];
            sum_k4 += (sGC.d.valbpt[0][k] + sGC.gamma_pt[k/2][sGC.time] + sGC.d.bpt[0][k]*sGC.delta_t[sGC.time])*v_k[k];
        }

        IloExpr obj(env);
        obj = sum_i2 - sum_k4;

        model.add(sum_k1 - sum_i1 <= sGC.p.qmax);
        model.add(sum_i1 - sum_k2 <= sGC.p.qmax);
        model.add(sum_k3 == 1);
        model.add(IloMaximize(env, obj));
        
		cplex.exportModel("cplexmodele.lp");
		
        if (!cplex.solve()){
            cout << "Failed to optimize LP." << endl;
            statusint = -1;
            cout << cplex.getStatus() << endl;
            cout << cplex.getCplexStatus() << endl;
            cout <<	cplex.getCplexSubStatus() << endl;
        }
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
            cplex.getValues(valUi,u_i);
            cplex.getValues(valVk,v_k);
            //if(pbdata.Params.CHOIX_AFFICHAGE==1){
                //cout << "Solution status = " << cplex.getStatus() << endl;
                //cout << "Solution value = " << cplex.getObjValue() << endl;
                //cout << "Values alpha = " << valX << endl;
            //}
        }
        if(statusint==0){
            if (SCIPisGT(sGC.scip,objfctvalue,0))
            {
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
                cout << "cout reduit SP CPLEX : " << objfctvalue << endl;
                status1 = addObjectColumnInModel(sGC,valUi,valVk,sGC.time);
                assert(status1==1);
                sGC.time++;
                //assert(pbdata.time - inittime >= 0);
                //cpteur = cpteur + (pbdata.time - inittime);
                //if(pbdata.Params.MULTIPLE_ENSEMBLE==0)
                stop = true;
            }else
            {   
                cout << "Cout reduit non négatif au temps "<<sGC.time<< endl;
                status = SCIP_SUCCESS;
                sGC.time++;
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
    // ??????????????? START SECTION TO DELETE
    /*if((!colgen)&&(sGC.todelete == 0)){
        IloEnv env;

        //ajout y320
        IloNumArray valUi3(env, 4, 0., 0., 1.,0.);
        IloNumArray valVk2(env, 4, 0., 0., 1.,0.);
        cout << "------------------ ADD COLONNE "<<endl;
        addObjectColumnInModel(sGC,valUi3,valVk2,0);
        
        //ajout y421
        IloNumArray valUi4(env, 4, 1., 1., 0.,0.);
        cout << "------------------ ADD COLONNE "<<endl;
        addObjectColumnInModel(sGC,valUi4,valVk2,1);
        
        //ajout y322 y332
        IloNumArray valVk3(env, 4, 0., 0., 0.,1.);
        cout << "------------------ ADD COLONNE "<<endl;
        addObjectColumnInModel(sGC,valUi3,valVk3,2);

        sGC.todelete = 1;
        env.end();
    }*/
    
 




    // ???????????????? END

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
	//pbdata->tempsPricer -= (double) clock()/CLOCKS_PER_SEC;
	//int duration = pbdata->instance.getDeadLine()-pbdata->instance.getReleaseTime();

	/*if(pbdata->Params.ACTIVATE_VERIF==1){
		assert(Pr_VerifSol(*pbdata,duration)==true);
	}*/
	//pbdata->attemptPrice = pbdata->attemptPrice + 1;
	// --------------------------------------   recuperation des variables duales
    int fois = 0;
    for(int i=0; i<pbdata->d.cardJ; ++i){
        pbdata->beta_i[i] = SCIPgetDualsolLinear(scip,pbdata->cons_2[i]);
		for(int t=0; t<pbdata->d.cardT; ++t){
            if((pbdata->d.rj[i]<=t)&&(t<pbdata->d.dj[i])){
                pbdata->alpha_it[i][t-pbdata->d.rj[i]] = SCIPgetDualsolLinear(scip,pbdata->cons_1[i][t-pbdata->d.rj[i]]);
            }
            if(fois==0){
                for(int p=0; p<(pbdata->d.nb_bp[0]/2); ++p){
                    pbdata->gamma_pt[p][t] = SCIPgetDualsolLinear(scip,pbdata->cons_3[p][t]);
                }
                if(t<pbdata->d.cardT-1) pbdata->delta_t[t] = SCIPgetDualsolLinear(scip,pbdata->cons_8[t]);
                pbdata->phi_t[t] = SCIPgetDualsolLinear(scip,pbdata->cons_9[t]);
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
	assert(status==1);
	pbdata->tempsPricer += (double) clock()/CLOCKS_PER_SEC;*/
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
        for(const auto& k : pbdata->K_l[l]){
		    for(int t=pbdata->L[l].releaseTime; t<pbdata->L[l].deadLine; ++t){
	            SCIPgetTransformedVar(scip, pbdata->varY_lkt[l][k][t],&(pbdata->varY_lkt[l][k][t]));
			    //assert(pbdata->varY_lkt[l][k][t] != NULL);
            }
		}
	}

    for(int k=0; k<pbdata->d.nb_bp[0]; ++k){
		for(int t=0; t<pbdata->d.cardT; ++t){
            SCIPgetTransformedVar(scip, pbdata->varY0_kt[k][t],&(pbdata->varY0_kt[k][t]));
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
                for(int p=0; p<(pbdata->d.nb_bp[0]/2); ++p){
                    SCIPgetTransformedCons(scip,pbdata->cons_3[p][t],&(pbdata->cons_3[p][t]));
                    SCIPgetTransformedCons(scip,pbdata->cons_4[p][t],&(pbdata->cons_4[p][t]));
                }
                //if(t<pbdata->d.cardT-1) 
                SCIPgetTransformedCons(scip,pbdata->cons_8[t],&(pbdata->cons_8[t]));
                SCIPgetTransformedCons(scip,pbdata->cons_9[t],&(pbdata->cons_9[t]));
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