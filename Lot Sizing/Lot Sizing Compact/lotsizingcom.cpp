#include <iostream>
#include <vector>
#include <string>
#include "scip/scip.h"
#include <scip/scipdefplugins.h>
using namespace std;

int main(){

	SCIP * scip;
	SCIP_CALL_EXC(SCIPcreate(& scip));
	
	return 0;
}

