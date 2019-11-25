#include<cassert>
#include<fstream>
#include<tuple>
#include<vector>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
using namespace std ;

int main( int argc, char* argv[]){    

    hopping_list hl_t;
    hopping_list::cellID_t cellID;
    hopping_list::edge_t  vertex_edge;
    hopping_list::value_t  hop;
    {
        //GRAPHENE TEST
        hopping_list hl_i = create_hopping_list(read_wannier_file("graphene_hr.dat"));
        hl_t = hopping_list();
        hl_t.num_wann = 2;
        cellID={ 0, 0, 0}; vertex_edge={0,1}; hop = complex<double>(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={0,1}; hop = complex<double>(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={0,1}; hop = complex<double>(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 0, 0}; vertex_edge={1,0}; hop = complex<double>(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={1,0}; hop = complex<double>(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={1,0}; hop = complex<double>(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={0,0}; hop = 1.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={0,0}; hop = 2.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={0,0}; hop = 3.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={0,0}; hop = complex<double>(4.123456,0.100000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 1, 0}; vertex_edge={0,0}; hop = complex<double>(5,-0.200000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1,-1, 0}; vertex_edge={0,0}; hop = complex<double>(1.0,1.0);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={1,1}; hop = 1.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={1,1}; hop = 2.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={1,1}; hop = 3.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={1,1}; hop = complex<double>(4.123456,-0.100000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 1, 0}; vertex_edge={1,1}; hop = complex<double>(4.123456,-0.100000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1,-1, 0}; vertex_edge={1,1}; hop = complex<double>(1.0,-1.0); 
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        assert( hl_t == hl_i );
    }
    //UNIT_CELL TEST   
    {
        hopping_list hl_i = create_hopping_list(read_wannier_file("graphene_hr.dat"));
        const array<int, 3> cellDim={1,1,1};
        hl_i = wrap_in_supercell(cellDim, hl_i );
        hl_t = hopping_list();
        hl_t.num_wann = 2;
        hl_t.cellSizes = cellDim;
        cellID={0,0,0}; vertex_edge={0,1}; hop = 3.0*complex<double>(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,0}; hop = 3.0*complex<double>(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={0,0}; hop = 1.0+2.0+3.0+ complex<double>(5,-0.200000) + complex<double>(4.123456,0.100000) + complex<double>(1.0,1.0);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,1}; hop = 1.0+2.0+3.0+ complex<double>(1.0,-1.0)+ complex<double>(4.123456,-0.100000) + complex<double>(4.123456,-0.100000);;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        assert( hl_t == hl_i );
    }
    /// SUPER_CELL TEST: Insert the hopping elements of a one dimensional chain 
    /// and construct a supercell using wrap_in_supercell function.
    /// compare this with a supercell constructed inserting the hoppings directly
    /// 
    {
        hopping_list 
        hl_t  =hopping_list(); //input
        cellID={ 1,0,0 }; vertex_edge={0,0}; hop = complex<double>(1.0, 1.0);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1,0,0 }; vertex_edge={0,0}; hop = complex<double>(1.0,-1.0);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        const array<int, 3> cellDim={2,1,1};


        hopping_list 
        hl_t = hopping_list();//test
        hl_t.num_wann = 2;
        hl_t.cellSizes = cellDim;
        cellID={0,0,0}; vertex_edge={0,1}; hop = 3.0*complex<double>(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,0}; hop = 3.0*complex<double>(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={0,0}; hop = 1.0+2.0+3.0+ complex<double>(5,-0.200000) + complex<double>(4.123456,0.100000) + complex<double>(1.0,1.0);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,1}; hop = 1.0+2.0+3.0+ complex<double>(1.0,-1.0)+ complex<double>(4.123456,-0.100000) + complex<double>(4.123456,-0.100000);;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        assert( hl_t == hl_i );
    }
 

return 0;}
