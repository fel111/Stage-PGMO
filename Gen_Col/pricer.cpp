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
#include "../struct.h"
#include "struct_gencol.h"
#include "pricer.h"

using namespace std;

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


/*! \brief Method creating and solving the subproblem for a predefined t
 It gives the subset of tasks in execution at time t (use a piece-wise linear obj function for the SP)
Then the column (subset of tasks) is added to the restricted master problem with one of the three adding methods
 *
 \param[in,out] pbdata Reference to structGenCol the entire problem data structure
 *
 * \return status: a variable of type SCIP_RESULT, equal to SCIP_SUCCESS if everything went well
 */
SCIP_RESULT Pr_SP1(structGenCol &sGC){
	
    double objfctvalue = 0;
    
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);

    //ajout variables u_i
	IloNumVarArray u_i (env,sGC.d.cardJ,0.0,1.0,ILOBOOL);
    
    //ajout variables v_k
	IloNumVarArray v_k (env,sGC.d.nb_bp[0],0.0,1.0,ILOBOOL);
    

    //contrainte sum_i bi*ui - sum_k o_k2(pk+1)*vk <= qmax
    IloExpr sum_i1(env);
    IloExpr sum_i2(env);
    for(int i=0; i<sGC.d.cardJ;++i){
        if((sGC.d.rj[i]<=sGC.time)&&(sGC.time<sGC.d.dj[i])){
            sum_i1 += sGC.d.Djk[i][0]*u_i[i];
            sum_i2 += (sGC.beta_i[i] - sGC.alpha_it[i][sGC.time] + (sGC.delta_t[sGC.time]+sGC.phi_t[sGC.time])*sGC.d.Djk[i][0])*u_i[i];
        }
        else model.add(u_i[i] <= 0);
    }
    IloExpr sum_k1(env);
    for(int k=0; k<sGC.d.nb_bp[0]-2; ++k){
        if(k%2==0) sum_k1 += sGC.d.bpt[0][k+3]*v_k[k];
        else sum_k1 += sGC.d.bpt[0][k+2]*v_k[k];
    }
    IloExpr sum_k2(env);
    for(int k=2; k<sGC.d.nb_bp[0]; ++k){
        if(k%2==0) sum_k2 += sGC.d.bpt[0][k-1]*v_k[k];
        else sum_k2 += sGC.d.bpt[0][k-2]*v_k[k];
    }
    IloExpr sum_k3(env);
    IloExpr sum_k4(env);
    for(int k=0; k<sGC.d.nb_bp[0]; ++k){
        sum_k3 += v_k[k];
        sum_k4 += (sGC.d.valbpt[0][k] + sGC.gamma_pt[k/2][sGC.time] + sGC.d.bpt[0][k]*sGC.delta_t[sGC.time])*v_k[k];
    }

    IloExpr obj(env);
    obj = sum_i2 - sum_k4;

    model.add(sum_i - sum_k1 <= sGC.p.qmax);
    model.add(sum_i - sum_k2 >= sGC.p.qmax);
    model.add(sum_k3 == 1);
    model.add(IloMaximize(env, obj));
    
    if (!cplex.solve()){
        cout << "Failed to optimize LP." << endl;
        cout << cplex.getStatus() << endl;
        cout << cplex.getCplexStatus() << endl;
        cout <<	cplex.getCplexSubStatus() << endl;
    }
    if(cplex.getStatus() !=IloAlgorithm::Optimal){
        cerr << "Cplex did not find the optimal solution for the subproblem !!!!!!!!!! " << endl;
        cout << "elapsed time=" << SCIPgetSolvingTime(sGC.scip) << endl;
        cout << cplex.getStatus() << endl;
        cout << cplex.getCplexStatus() << endl;
        cout <<	cplex.getCplexSubStatus() << endl;
        cout << "obj =" << cplex.getObjValue() << endl;
        statusint = -1;
    }
    }else{
        objfctvalue = cplex.getObjValue();
        cplex.getValues(valX,varAlpha);
        //if(pbdata.Params.CHOIX_AFFICHAGE==1){
            cout << "Solution status = " << cplex.getStatus() << endl;
            cout << "Solution value = " << cplex.getObjValue() << endl;
            cout << "Values alpha = " << valX << endl;
        //}
    }
    
    
    
    
    
    
    
    
    
    // -------------------------
    // old code
    
    bool forcexit = false;
	int statusint = 0;
	int cpteur = 0;
	int inittime = 0;
	double objfctvalue = 0;
    bool colgen=false;

	SCIP_RESULT status = SCIP_SUCCESS;
	SCIP_RETCODE status1 = SCIP_ERROR;
	if(pbdata.Params.INCREMENTATION_TEMPS==0)
		pbdata.time = 0;
	else
		if(pbdata.time == pbdata.instance.getDeadLine())
			pbdata.time = 0;
	while(cpteur <= pbdata.instance.getDeadLine() && forcexit == false){
		
		
			
			
			pbdata.tempsResol += (double) clock()/CLOCKS_PER_SEC;
			
			
		}

		catch (IloException& e) {
			cerr << "Concert exception caught: " << e << endl;
			statusint=-1;
			colgen=false;
		}
		catch (...) {
			cerr << "Unknown exception caught" << endl;
			statusint=-1;
			colgen=false;
		}
		if(statusint==0){
			if (SCIPisGT(pbdata.scip,objfctvalue,0.0))
			{
				inittime = pbdata.time;
				// set result pointer
                                colgen=true;
				status = SCIP_SUCCESS;
				switch (pbdata.Params.TYPE_AJOUT)
				{
				case 1:
					status1 = Pr_addObjectColumnInModel_1(pbdata,valX,valY,valZ,objfctvalue);
					break;
				case 2:
					status1 = Pr_addObjectColumnInModel_2(pbdata,valX,valY,valZ,objfctvalue);
					break;
				case 3:
					status1 = Pr_addObjectColumnInModel_3(pbdata,valX,valY,valZ,objfctvalue);
					break;
				}
				assert(status1==1);
				assert(pbdata.time - inittime >= 0);
				cpteur = cpteur + (pbdata.time - inittime);
				if(pbdata.Params.MULTIPLE_ENSEMBLE==0)
					forcexit = true;
			}else
			{
				status = SCIP_SUCCESS;
				pbdata.time++;
				cpteur = cpteur + 1;
			}
		}else{
			cout << "colgen=" << colgen << endl;
                        if (!colgen)
                                status = SCIP_DIDNOTRUN;
                        else
                                status=SCIP_SUCCESS;

			forcexit = true;
		}
		if(pbdata.time == pbdata.instance.getDeadLine())
			pbdata.time = 0;
		env.end();
	}
	//cout << "end SP" << endl;
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
            if((t<=pbdata->d.rj[i])&&(t<pbdata->d.dj[i])){
                pbdata->alpha_it[i][t] = SCIPgetDualsolLinear(scip,pbdata->cons_1[i][t]);
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
	            SCIPgetTransformedVar(scip, pbdata->varY_lkt[l][k][t],&(pbdata->varY_lkt[l][k][t]);
			    //assert(pbdata->varY_lkt[l][k][t] != NULL);
            }
		}
	}

	for(int i=0; i<pbdata->d.cardJ; i++){
		for(int t = pbdata->d.rj[i] ; t < pbdata->d.dj[i] ; t++){
			SCIPgetTransformedVar(scip, pbdata->varX_it[i][t],&varX_it[i][t]);
			//assert(pbdata->tabVarXit.at(i).at(j-pbdata->instance.getTasksList().at(i).getReleaseTime()) != NULL);
		}
	}

	// get transformed constraints
	//int duration = pbdata->instance.getDeadLine()-pbdata->instance.getReleaseTime();
	int fois = 0;

	for(int i=0; i<pbdata->d.cardJ; ++i){
        SCIPgetTransformedCons(scip,pbdata->cons_2[i],&(pbdata->cons_2[i]));
        
		for(int t=0; t<pbdata->d.cardT; ++t){
            if((t<=pbdata->d.rj[i])&&(t<pbdata->d.dj[i])){
                SCIPgetTransformedCons(scip,pbdata->cons_1[i][t],&(pbdata->cons_1[i][t]));
                //assert(pbdata->tabConsNbSets[j] != NULL);
            }
            if(fois==0){
                for(int p=0; p<(pbdata->d.nb_bp[0]/2); ++p){
                    SCIPgetTransformedCons(scip,pbdata->cons_3[p][t],&(pbdata->cons_3[p][t]));
                }
                if(t<pbdata->d.cardT-1) SCIPgetTransformedCons(scip,pbdata->cons_8[t],&(pbdata->cons_8[t]));
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