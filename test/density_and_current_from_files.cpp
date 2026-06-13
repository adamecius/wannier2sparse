#include<cassert>
#include<fstream>
#include<tuple>
#include<vector>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
#include "tbmodel.hpp"
using namespace std ;

// Read-only "value at (cellID, edge)" lookup, summing any matching terms and
// returning zero when the key is absent.  The flat hopping_storage::operator[]
// asserts on a missing key, so this restores the map-style read semantics these
// tests rely on (a missing contribution counts as zero) without mutating hl.
static hopping_list::value_t value_at(const hopping_list& hl,
                                      const hopping_list::cellID_t& cid,
                                      const hopping_list::edge_t& e)
{
    hopping_list::value_t sum(0.0, 0.0);
    for (const auto& h : hl.hoppings)
        if (get<0>(h) == cid && get<2>(h) == e)
            sum += get<1>(h);
    return sum;
}

int main( int argc, char* argv[]){    

    tbmodel model;

    model.readOrbitalPositions("spin_graphene.xyz");
    model.readUnitCell("spin_graphene.uc");
    model.readWannierModel("spin_graphene_hr.dat");    

    hopping_list JxSx= model.createHoppingSpinCurrents_list(0,'x');
    hopping_list JxSy= model.createHoppingSpinCurrents_list(0,'y');
    hopping_list JxSz= model.createHoppingSpinCurrents_list(0,'z');

    hopping_list Sx = model.createHoppingSpinDensity_list('x');
    hopping_list Sy = model.createHoppingSpinDensity_list('y');
    hopping_list Sz = model.createHoppingSpinDensity_list('z');

    hopping_list hl_i = model.hl;
    hopping_list hl_t;
    hopping_list::cellID_t cellID;
    hopping_list::edge_t  vertex_edge;
    hopping_list::value_t  hop;
    string tag ; 

    std::cout<<"Performing  WRAP IN SUPERCELL FOR UNIT_CELL TEST"<<std::endl;
    {
        const hopping_list::cellID_t SCDIM({1,1,1});
        const int unit_cell_basis_size = 4;

        hl_t = hopping_list();
        hl_t.SetWannierBasisSize(unit_cell_basis_size);
        hl_t.cellSizes = SCDIM;


        cellID={0 , 0, 0}; vertex_edge={0,2};
        hop  = value_at(hl_i, cellID, vertex_edge);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={2,0};
        hop  = value_at(hl_i, cellID, vertex_edge);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );


        cellID={0 , 0, 0}; vertex_edge={0,1};
        hop  = value_at(hl_i, cellID, vertex_edge);
        cellID={-1, 0, 0}; vertex_edge={0,1};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={ 0,-1, 0}; vertex_edge={0,1};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={1,0};
        hop  = value_at(hl_i, cellID, vertex_edge);
        cellID={ 1, 0, 0}; vertex_edge={1,0};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={ 0, 1, 0}; vertex_edge={1,0};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={2,3};
        hop  = value_at(hl_i, cellID, vertex_edge);
        cellID={-1, 0, 0}; vertex_edge={2,3};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={ 0,-1, 0}; vertex_edge={2,3};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={3,2};
        hop  = value_at(hl_i, cellID, vertex_edge);
        cellID={ 1, 0, 0}; vertex_edge={3,2};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={ 0, 1, 0}; vertex_edge={3,2};
        hop += value_at(hl_i, cellID, vertex_edge);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );


        // The diagonal on-site terms are exactly zero in spin_graphene_hr.dat,
        // so create_hopping_list drops them and the wrapped supercell must not
        // contain them.  Assert their absence instead of inserting zero entries
        // (the old map operator[] used to fabricate matching zeros as a side
        // effect of the lookup; the flat storage does not).
        assert( value_at(hl_i, hopping_list::cellID_t({0,0,0}), hopping_list::edge_t({0,0})) == hopping_list::value_t(0.0,0.0) );
        assert( value_at(hl_i, hopping_list::cellID_t({0,0,0}), hopping_list::edge_t({1,1})) == hopping_list::value_t(0.0,0.0) );
        assert( value_at(hl_i, hopping_list::cellID_t({0,0,0}), hopping_list::edge_t({2,2})) == hopping_list::value_t(0.0,0.0) );
        assert( value_at(hl_i, hopping_list::cellID_t({0,0,0}), hopping_list::edge_t({3,3})) == hopping_list::value_t(0.0,0.0) );


        hl_i = wrap_in_supercell(SCDIM, hl_i );
        assert( hl_t == hl_i );
    }
    std::cout<<"WRAP IN SUPERCELL FOR UNIT_CELL TEST"<<std::endl<<std::endl;

return 0;}
