#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <deque>
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

hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  )
{
    const int num_wann= get<0>(wannier_data);
    const vector<string> hopping_lines= get<1>(wannier_data);
    vector< tuple<int,int> > vertex_edge;
    map<string, tuple< array<int, 3>,complex<double> > >hoppings; 


    stringstream ss;
    for (auto line : hopping_lines){
        ss.str(line);
        
        int forb,torb; 
        ss>>forb>>torb; 
        vertex_edge.push_back(tuple<int,int>(forb,torb) );
        
        array<int, 3> tag; string stag; 
        for( auto& ti : tag){
            ss>>ti;
            stag += to_string(ti)+" ";
        }
        double rv,im; ss>>rv>>im;
        hoppings.insert( {stag,tuple(tag,complex(rv,im) ) } );
    }
    
hopping_list output;
output.num_wann =num_wann;
output.vertex_edge = vertex_edge;
output.hoppings = hoppings;

return output; }



int main( int argc, char* argv[]){

deque< string > arguments(argv,argv+argc);
const string program_name = arguments[0]; arguments.pop_front();  

if( arguments.empty() )
    std::cerr<<"ERROR: The program: "<<program_name <<" should be called with at least one argument (LABEL). "<<std::endl;

const string  label = arguments[0]; arguments.pop_front();  

std::cout<<"Using "<<system<<" as the system's identification label"<<std::endl
         <<"This label will be used to detect the label.xyz, label_hr.dat, and label.win files"<<std::endl;

const string wannier_filename = label+"_hr.dat";

hopping_list hops  = create_hopping_list(read_wannier_file(wannier_filename));


return 0;}
