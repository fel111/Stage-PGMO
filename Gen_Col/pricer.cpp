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


double verifCoutReduit(structGenCol const& sGC, vector<int> ui, int k, int t){
    double cr = 0.0;
    for(int i=0; i<ui.size(); ++i){
        if(SCIPisEQ(sGC.scip, ui[i], 1.0)) cr += (sGC.beta_i[i] - sGC.alpha_it[i][t-sGC.d.rj[i]] + sGC.d.Dj[i]*(sGC.delta_t[t] + sGC.phi_t[t]));
    }
    cr -= (sGC.d.valbpt[0][k] + sGC.gamma_pt[k/2][t] + sGC.d.bpt[0][k]*sGC.delta_t[t]);
    return -cr;
}

double verifCoutReduit(structGenCol const& sGC, IloNumArray ui, int k, int t){
    double cr = 0.0;
    for(int i=0; i<ui.getSize(); ++i){
        if(SCIPisEQ(sGC.scip, ui[i], 1.0)) cr += (sGC.beta_i[i] - sGC.alpha_it[i][t-sGC.d.rj[i]] + sGC.d.Dj[i]*(sGC.delta_t[t] + sGC.phi_t[t]));
    }
    cr -= (sGC.d.valbpt[0][k] + sGC.gamma_pt[k/2][t] + sGC.d.bpt[0][k]*sGC.delta_t[t]);
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
SCIP_RETCODE addObjectColumnInModel (structGenCol &sGC,IloNumArray valUi,IloNumArray valVk,int timedd){
	SCIP_RETCODE status = SCIP_OKAY; // devrait etre modifie lors de l ajout d une colonne, sinon erreur
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
    
    int k = 0 ;
    while( SCIPisEQ(sGC.scip, valVk[k], 0.0) ) ++k;


    if(taskList.size() > 0){
        l.id = sGC.L.size();
        l.timeGen = timedd; 
        l.tasksList = taskList;
        int testPres = checkSet(l,sGC);
        int id;
        bool createOk = false;
        if(testPres == -1){
            sGC.L.push_back(l);
            sGC.cardL = sGC.L.size();
            addSetK_l(l,sGC);
            addA_il(l,sGC);
            addL_t(l,sGC);
            id = sGC.L.size()-1;
            createOk = true;
        }
        else{
            id = testPres;
            if(sGC.varY_lkt[id][k][timedd] == NULL) createOk = true;
        }
        if(createOk){
            // Creation de la nouvelle variable pour le moment t obtenu a l issu du sous probleme
            string name = "y_lkt"+to_string(id)+","+to_string(k)+","+to_string(timedd);
            //cout <<"ajout variable "+name<<endl;
            SCIP_VAR * var;
            SCIPcreateVarBasic(sGC.scip, &var, name.c_str() ,0,SCIPinfinity(sGC.scip),sGC.d.valbpt[0][k],SCIP_VARTYPE_CONTINUOUS);
            //SCIPcreateVarBasic(sGC.scip, &var, NULL ,0,SCIPinfinity(sGC.scip),sGC.d.valbpt[0][k],SCIP_VARTYPE_CONTINUOUS);
            for(const auto& j : taskList){
                //Mise a jour cons1
                //cout<<"j : "<<j<<" timedd : "<<sGC.timedd<<" rj : "<<sGC.d.rj[j] << " dj : "endl;
                //SCIPprintCons(sGC.scip,sGC.cons_1[j][sGC.time-sGC.d.rj[j]],NULL);
                SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_1[j][timedd-sGC.d.rj[j]],var,-1));
                // Mise a jour cons2
                SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_2[j],var,1));
            }
            // Mise a jour cons3
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_3[k/2][timedd],var,-1));
            // Mise a jour cons8
            //cout<<"contrainte 8 : l.energyDemand-sGC.d.bpt[0][k] = "<<l.energyDemand-sGC.d.bpt[0][k]<<endl;
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_8[timedd],var,l.energyDemand-sGC.d.bpt[0][k]));
            // Mise a jour cons9
            SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_9[timedd],var,l.energyDemand));
            SCIP_CALL(status = SCIPaddPricedVar(sGC.scip, var, 1.0));
            //if(sGC.todelete2 == 0) sGC.todelete3 = SCIPgetVarRedcost(sGC.scip,var);
            //assert( SCIPisLT(sGC.scip,SCIPgetVarRedcost(sGC.scip,var),0) );
            //cout << "valeur cout reduit -a la main- : " << verifCoutReduit(sGC,valUi,k,timedd) << "  valeur cout reduit scip : " << SCIPgetVarRedcost(sGC.scip,var) << endl;
            //assert(SCIPisEQ(sGC.scip,verifCoutReduit(sGC,valUi,k,timedd),SCIPgetVarRedcost(sGC.scip,var)));
            if(testPres==-1){
                vector<vector<SCIP_VAR*> > tab (sGC.d.nb_bp[0], vector<SCIP_VAR*> (sGC.d.cardT, NULL));
                tab[k][timedd] = var;
                sGC.varY_lkt.push_back(tab);
            }
            else{
                sGC.varY_lkt[id][k][timedd] = var;
            }
            sGC.nbcolgenerated++;
        
            //else cout << "colonne existante !!!!!!!"  << endl;
            // PARTIE AJOUT COLONNE Y_LHT, AVEC H LE BREAKPOINT VOISIN DE K
            
            if(sGC.p.ajout_breakpoint_voisin == 1){
                int k2;
                if(k%2==0) k2=k+1;
                else k2=k-1; 
                if((testPres==-1)||((testPres!=-1)&&(sGC.varY_lkt[id][k2][timedd]==NULL))){
                    string name2 = "y_lkt"+to_string(id)+","+to_string(k2)+","+to_string(timedd);
                    SCIP_VAR * var2;
                    SCIPcreateVarBasic(sGC.scip, &var2, name2.c_str() ,0,SCIPinfinity(sGC.scip),sGC.d.valbpt[0][k2],SCIP_VARTYPE_CONTINUOUS);
                    for(const auto& j : taskList){
                        //Mise a jour cons1
                        //cout<<"j : "<<j<<" time : "<<sGC.time<<" rj : "<<sGC.d.rj[j] << " dj : "endl;
                        //SCIPprintCons(sGC.scip,sGC.cons_1[j][sGC.time-sGC.d.rj[j]],NULL);
                        SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_1[j][timedd-sGC.d.rj[j]],var2,-1));
                        // Mise a jour cons2
                        SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_2[j],var2,1));
                    }
                    // Mise a jour cons3
                    SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_3[k2/2][timedd],var2,-1));
                    // Mise a jour cons8
                    //cout<<"contrainte 8 : l.energyDemand-sGC.d.bpt[0][k] = "<<l.energyDemand-sGC.d.bpt[0][k]<<endl;
                    SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_8[timedd],var2,l.energyDemand-sGC.d.bpt[0][k2]));
                    // Mise a jour cons9
                    SCIP_CALL(status = SCIPaddCoefLinear(sGC.scip,sGC.cons_9[timedd],var2,l.energyDemand));
                    SCIP_CALL(status = SCIPaddPricedVar(sGC.scip, var2, 1.0));
                    sGC.nbcolgenerated++;
                    // ???????????? START
                    //cout << "cout reduit de colonne ajoutee : "<<SCIPgetVarRedcost(sGC.scip,var2)<<endl;
                }
            }
            //return status;
        }
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
    return status;
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
	
    //cout << "dual bound : " << SCIPgetDualbound(sGC.scip) << endl;
    //cout << "nodes left : " << SCIPgetNNodesLeft(sGC.scip) << endl;
    //cout << "nodes : " << SCIPgetNNodes(sGC.scip) << endl;

    /*if(SCIPgetDualbound(sGC.scip)>=229.8) sGC.todelete++;
    if(sGC.todelete == 8){
        FILE * filed;
        filed = fopen("before230_85","w");
        SCIPprintTransProblem(sGC.scip, filed, "lp", 0);
        fclose(filed);
        /*
        for(int l=0; l<sGC.L.size(); ++l){
            for(const auto& k : sGC.K_l[l]){
                for(int t=sGC.L[l].releaseTime; t<sGC.L[l].deadLine; ++t){
                    if(sGC.varY_lkt[l][k][t] != NULL){
                        //if(SCIPisGT(sGC.scip,SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY_lkt[l][k][t]),0) ) cout << "y"<<l<<","<<k<<","<<t<<" = "<<SCIPgetSolVal(sGC.scip, sGC.sol,sGC.varY_lkt[l][k][t])<<endl;
                        cout << "y"<<l<<","<<k<<","<<t<<endl;
                    }
                }
            }
        }* /
	
    }*/

    int timedd = 0;
    double objfctvalue = 0.0;
    int statusint = 0;
    bool colgen = false;
    bool stop = false;
    SCIP_RESULT status = SCIP_SUCCESS;
	SCIP_RETCODE status1 = SCIP_ERROR;
    
    while((timedd<sGC.d.cardT)&&(!stop)){

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
			//cout << i << "    " << sGC.d.rj[i] << " <= " << timedd<<" && "<< timedd<<" < " <<sGC.d.dj[i]<<endl;
            if((sGC.d.rj[i]<=timedd)&&(timedd<sGC.d.dj[i])){
                sum_i1 += sGC.d.Dj[i]*u_i[i]; //sum bi ui
                //???
                //cout << "alpha : "<<sGC.alpha_it[i] << endl;
                //cout << "time" << timedd << endl;
                sum_i2 += u_i[i] * (sGC.beta_i[i] //sum_i obj
                    - sGC.alpha_it[i][timedd-sGC.d.rj[i]] 
                    + sGC.d.Dj[i] 
                        * (sGC.delta_t[timedd]
                           + sGC.phi_t[timedd]));
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
            sum_k4 += (sGC.d.valbpt[0][k] + sGC.gamma_pt[k/2][timedd] + sGC.d.bpt[0][k]*sGC.delta_t[timedd])*v_k[k];
        }

        IloExpr obj(env);
        obj = sum_i2 - sum_k4;

        model.add(sum_k1 - sum_i1 <= sGC.d.Q);
        model.add(sum_i1 - sum_k2 <= sGC.d.Q);
        model.add(sum_k3 == 1);
        model.add(IloMaximize(env, obj));
        
		//cplex.exportModel("cplexmodele.lp");
		
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
                //cout << "cout reduit SP CPLEX : " << -objfctvalue << endl;
                status1 = addObjectColumnInModel(sGC,valUi,valVk,timedd);
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
    //if(!colgen) cout << "no colonnes !!" << endl;
    // ??????????????? START SECTION TO DELETE

    /*if(!colgen){
        IloEnv env;
        IloNumArray valUi1(env, 10, 1, 1, 0,1,0,1,0,1,0,0);
        cout << "cout reduit colonne ameliorante : " << verifCoutReduit(sGC,valUi1,3,1);
    
        if((SCIPnodeGetNumber(SCIPgetCurrentNode(sGC.scip)) >= 58.9)&&(SCIPnodeGetNumber(SCIPgetCurrentNode(sGC.scip)) <= 59.1)){
            cout << "good node ! " << endl;
            SCIPwriteMIP(sGC.scip,"mipnode",0,1,1);
            SCIPwriteLP(sGC.scip,"lpnode");
            FILE *filed;
            filed = fopen("transprobnode","w");
            SCIPprintTransProblem(sGC.scip,filed,"lp",0);
            fclose(filed);
        }
    }*/
    
    //sGC.todelete = 0;

    /*if((!colgen)&&(sGC.todelete2 == 0)){
        //FILE * filed;
        //filed = fopen("beforeCol", "a+");
        //cout << " BEFORE COL ------------------------- " << endl;
        //SCIPprintTransProblem(sGC.scip, NULL, "lp", 0);
        //fclose(filed);

        IloEnv env;

        //ajout y320
        IloNumArray valUi1(env, 10, 1, 1, 0,1,0,1,0,1,0,0);
        IloNumArray valUi2(env, 10, 1,0,0,0,0,0,0,1,1,0);
        IloNumArray valUi3(env, 10, 1,0,0,0,1,1,0,0,1,1);
        IloNumArray valUi4(env, 10, 1,0,0,0,0,0,0,0,0,0);
        IloNumArray valUi5(env, 10, 0,0,1,0,0,1,1,0,0,0);
        IloNumArray valUi6(env, 10, 1,0,1,0,0,1,1,0,1,1);
        IloNumArray valUi7(env, 10, 1,0,1,0,0,0,0,0,0,0);
        IloNumArray valUi8(env, 10, 0,0,1,0,0,0,1,0,0,0);
        IloNumArray valUi9(env, 10, 0,0,1,0,0,0,0,0,0,0);
        IloNumArray valVk1(env, 4, 1., 0., 0.,0.);
        IloNumArray valVk2(env, 4, 0., 1., 0.,0.);
        IloNumArray valVk3(env, 4, 0., 0., 1.,0.);
        IloNumArray valVk4(env, 4, 0., 0., 0.,1.);
        cout << "------------------ ADD COLONNE "<<endl;
        //addObjectColumnInModel(sGC,valUi1,valVk3,1);
        //addObjectColumnInModel(sGC,valUi1,valVk4,1);
        /*addObjectColumnInModel(sGC,valUi2,valVk1,2);
        addObjectColumnInModel(sGC,valUi2,valVk2,2);
        addObjectColumnInModel(sGC,valUi3,valVk3,3);
        addObjectColumnInModel(sGC,valUi3,valVk4,3);
        /*addObjectColumnInModel(sGC,valUi4,valVk1,4);
        addObjectColumnInModel(sGC,valUi4,valVk2,4);
        addObjectColumnInModel(sGC,valUi5,valVk1,5);
        addObjectColumnInModel(sGC,valUi5,valVk2,5);
        addObjectColumnInModel(sGC,valUi6,valVk3,6);
        addObjectColumnInModel(sGC,valUi6,valVk4,6);
        addObjectColumnInModel(sGC,valUi7,valVk1,7);
        addObjectColumnInModel(sGC,valUi7,valVk2,7);
        addObjectColumnInModel(sGC,valUi8,valVk1,8);
        addObjectColumnInModel(sGC,valUi8,valVk2,8);
        addObjectColumnInModel(sGC,valUi9,valVk1,9);
        addObjectColumnInModel(sGC,valUi9,valVk2,9);* /


        sGC.todelete2 = 1;
        env.end();
        /*
        //cout << " AFTER COL ------------------------- " << endl;
        FILE * filed;
        filed = fopen("modeleAfterCol","w");
        SCIPprintTransProblem(sGC.scip, filed, NULL, 0);
        fclose(filed);* /
    }*/
    
 




    // ???????????????? END

	return status;
}

SCIP_RESULT Pr_farkas(structGenCol &sGC){
    int timedd = 0;
    double objfctvalue = 0.0;
    int statusint = 0;
    bool colgen = false;
    bool stop = false;
    SCIP_RESULT status = SCIP_SUCCESS;
	SCIP_RETCODE status1 = SCIP_ERROR;
    
    while((timedd<sGC.d.cardT)&&(!stop)){

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
			//cout << i << "    " << sGC.d.rj[i] << " <= " << timedd<<" && "<< timedd<<" < " <<sGC.d.dj[i]<<endl;
            if((sGC.d.rj[i]<=timedd)&&(timedd<sGC.d.dj[i])){
                sum_i1 += sGC.d.Dj[i]*u_i[i]; //sum bi ui
                //???
                //cout << "alpha : "<<sGC.alpha_it[i] << endl;
                //cout << "time" << timedd << endl;
                sum_i2 += u_i[i] * (sGC.farkas2_i[i] //sum_i obj
                    - sGC.farkas1_it[i][timedd-sGC.d.rj[i]] 
                    + sGC.d.Dj[i] 
                        * (sGC.farkas8_t[timedd]
                           + sGC.farkas9_t[timedd]));
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
            sum_k4 += (sGC.farkas3_pt[k/2][timedd] + sGC.d.bpt[0][k]*sGC.farkas8_t[timedd])*v_k[k];
        }

        IloExpr obj(env);
        obj = sum_i2 - sum_k4;

        model.add(sum_k1 - sum_i1 <= sGC.d.Q);
        model.add(sum_i1 - sum_k2 <= sGC.d.Q);
        model.add(sum_k3 == 1);
        model.add(IloMaximize(env, obj));
        
		//cplex.exportModel("cplexmodele.lp");
		
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
                //cout << "cout reduit SP CPLEX : " << -objfctvalue << endl;
                status1 = addObjectColumnInModel(sGC,valUi,valVk,timedd);
                assert(status1==1);
                timedd++;
                //assert(pbdata.time - inittime >= 0);
                //cpteur = cpteur + (pbdata.time - inittime);
                //if(pbdata.Params.MULTIPLE_ENSEMBLE==0)
                stop = true;
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
    //if(!colgen) cout << "no colonnes !!" << endl;


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
        pbdata->beta_i[i] = SCIPgetDualsolLinear(scip,pbdata->cons_2[i]);
		for(int t=0; t<pbdata->d.cardT; ++t){
            if((pbdata->d.rj[i]<=t)&&(t<pbdata->d.dj[i])){
                pbdata->alpha_it[i][t-pbdata->d.rj[i]] = SCIPgetDualsolLinear(scip,pbdata->cons_1[i][t-pbdata->d.rj[i]]);
            }
            if(fois==0){
                for(int p=0; p<(pbdata->d.nb_bp[0]/2); ++p){
                    pbdata->gamma_pt[p][t] = SCIPgetDualsolLinear(scip,pbdata->cons_3[p][t]);
                }
                //if(t<pbdata->d.cardT-1) pbdata->delta_t[t] = SCIPgetDualsolLinear(scip,pbdata->cons_8[t]);
                pbdata->delta_t[t] = SCIPgetDualsolLinear(scip,pbdata->cons_8[t]);
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
	assert(status==1);*/
	auto end_time = chrono::steady_clock::now();
    pbdata->tpsPricer += chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()/1000000.0;

	return SCIP_OKAY;
}


static
SCIP_DECL_PRICERFARKAS(pricerFarkas){
    structGenCol* pbdata; /* the data of the pricer */
	SCIP_PROBDATA* probdata;
	//SCIP_RETCODE status = SCIP_OKAY;
	assert(scip != NULL);
	assert(pricer != NULL);
	probdata = SCIPgetProbData(scip);
	pbdata = (structGenCol*) probdata;

    pbdata->cptFarkas++;
    //cout << "cptFarkas : " << pbdata->cptFarkas << endl;

    int fois = 0;
    for(int i=0; i<pbdata->d.cardJ; ++i){
        pbdata->farkas2_i[i] = SCIPgetDualfarkasLinear(scip,pbdata->cons_2[i]);
		for(int t=0; t<pbdata->d.cardT; ++t){
            if((pbdata->d.rj[i]<=t)&&(t<pbdata->d.dj[i])){
                pbdata->farkas1_it[i][t-pbdata->d.rj[i]] = SCIPgetDualfarkasLinear(scip,pbdata->cons_1[i][t-pbdata->d.rj[i]]);
            }
            if(fois==0){
                for(int p=0; p<(pbdata->d.nb_bp[0]/2); ++p){
                    pbdata->farkas3_pt[p][t] = SCIPgetDualfarkasLinear(scip,pbdata->cons_3[p][t]);
                }
                //if(t<pbdata->d.cardT-1) pbdata->delta_t[t] = SCIPgetDualfarkasLinear(scip,pbdata->cons_8[t]);
                pbdata->farkas8_t[t] = SCIPgetDualfarkasLinear(scip,pbdata->cons_8[t]);
                pbdata->farkas9_t[t] = SCIPgetDualfarkasLinear(scip,pbdata->cons_9[t]);
            }
        }
        ++fois;
	}

    *result = Pr_farkas(*pbdata);
		
	if (*result==SCIP_DIDNOTRUN) {
        cout << "DID NOT RUN " << endl;
        SCIPinterruptSolve(pbdata->scip);	
    }

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
                    //SCIPgetTransformedCons(scip,pbdata->cons_4[p][t],&(pbdata->cons_4[p][t]));
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
	SCIP_CALL(status =  SCIPincludePricerBasic(sGC.scip, &pricer, PRICER_NAMESP1, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY, pricerRedcost, pricerFarkas, pricerdata));
	assert(pricer != NULL);
	SCIP_CALL(status = SCIPsetPricerInitsol(sGC.scip, pricer, pricerInitsolSP));
    return status;
}