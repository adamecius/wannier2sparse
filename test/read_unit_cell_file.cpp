#include<cassert>
#include<fstream>
#include "wannier_parser.hpp"

int main( int argc, char* argv[]){

    array< array<double,3> , 3 >  
    unit_cell={
        array<double,3>({3.1416,1.1111,2.2222}),
        array<double,3>({1,2,3}),
        array<double,3>({0,1.,3e-2})
    };   
    const string unitcell_filename = "a_unit_cell.uc";

    ofstream unit_cell_file(unitcell_filename);
    for(const auto& lat_vec : unit_cell ){
        for(const auto& vec_elem : lat_vec )
            unit_cell_file<<vec_elem<<" ";
        unit_cell_file<<endl;
    }unit_cell_file.close();
    
    assert( unit_cell == read_unit_cell_file(unitcell_filename) );
    

return 0;}
