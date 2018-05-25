#include <vector>
#include <iostream>
#include "struct_gencol.h"
using namespace std;


float cl(int k, int l, structGenCol const& sGC){
    return sGC.d.valbpt[k];
}