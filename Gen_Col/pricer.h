/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#ifndef __SCIP_PRICER_H__
#define __SCIP_PRICER_H__

//#include "objscip/objscip.h"
//#include "scip/pub_var.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include <vector>
//#include <list>

using namespace std;
using namespace scip;


/** pricer class */
class ObjPricerPGMO : public ObjPricer
{
public:

   /** Constructs the pricer object with the data needed */
   ObjPricerPGMO(
      SCIP*                               scip,        /**< SCIP pointer */
      //const char*                         p_name,      /**< name of pricer */
      //const int                           p_num_nodes, /**< number of nodes */
      //const int                           p_capacity,  /**< vehicle capacity */
      //const vector< int >&                p_demand,    /**< demand array */
     // const vector< vector<int> >&        p_distance,  /**< matrix of distances */
     // const vector< vector<SCIP_VAR*> >&  p_arc_var,   /**< matrix of arc variables */
     // const vector< vector<SCIP_CONS*> >& p_arc_con,   /**< matrix of arc constraints */
     // const vector<SCIP_CONS* >&          p_part_con   /**< array of partitioning constraints */
	int nbVar;
	vector<SCIP_VAR*> var;
	vector<SCIP_CONS*> cons;
	vector<vector<int> > coeff;
	      

	);

   /** Destructs the pricer object. */
   virtual ~ObjPricerVRP();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** farkas pricing method of variable pricer for infeasible LPs */
   //virtual SCIP_DECL_PRICERFARKAS(scip_farkas);

};

#endif
