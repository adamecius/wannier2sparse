#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>

using namespace std;

struct hopping_list
{
    int num_wann;
    vector< tuple<int,int> > vertex_edge;
    map<string, tuple< array<int, 3>,complex<double> > > hoppings; 
};

hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  );
