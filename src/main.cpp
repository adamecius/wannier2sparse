#include <deque>
#include <string>
#include <iostream>

#include "wannier_parser.hpp"
#include "hopping_list.hpp"

using namespace std;

int main( int argc, char* argv[]){

deque< string > arguments(argv,argv+argc);
const string program_name = arguments[0]; arguments.pop_front();  

if( arguments.empty() )
{
    cerr<<"ERROR: The program: "<<program_name <<" should be called with at least one argument (LABEL). "<<endl;
    return -1;
}

const string  label = arguments[0]; arguments.pop_front();  

cout<<"Using "<<label<<" as the system's identification label"<<endl
    <<"This label will be used to detect the label.xyz, label_hr.dat, and label.win files"<<endl;

const string xyz_filename = label+".xyz";
read_xyz_file(xyz_filename);

const string uc_filename = label+".uc";
read_unit_cell_file(uc_filename);

const string wannier_filename = label+"_hr.dat";
const array<int, 3> cellDim={3,3,3};
hopping_list hl  =wrap_in_supercell(cellDim, create_hopping_list(read_wannier_file(wannier_filename) ) );

save_hopping_list_as_csr("test", hl);
cout<<"The programa finished"<<std::endl;
return 0;}
