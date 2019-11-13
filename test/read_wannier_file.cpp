#include<limits>
#include<iomanip>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<tuple>
#include<vector>

#include "wannier_parser.hpp"
using namespace std ;

int main( int argc, char* argv[]){    


    //GRAPHENE TEST
    int num_wann = 2;   
    vector<string> hopping_data;
    hopping_data.push_back("   0    0    0    1    1    0.000000   -0.000000");
    hopping_data.push_back("   0    0    0    2    2    0.000000   -0.000000");
    hopping_data.push_back("  -1    0    0    1    2    1.000000   -0.000000");
    hopping_data.push_back("   0   -1    0    1    2    1.000000   -0.000000");
    hopping_data.push_back("   1    0    0    2    1    1.000000   -0.000000");
    hopping_data.push_back("   0    1    0    2    1    1.000000   -0.000000");
    tuple<int, vector<string> > wan_data = read_wannier_file("graphene_hr.dat");

    assert(
        num_wann == get<0>(wan_data) &&
        hopping_data == get<1>(wan_data) 
        );

return 0;}
