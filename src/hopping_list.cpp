#include "hopping_list.hpp"
#include <functional>
#include<iostream>
#include<limits>


hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  )
{
    hopping_list hl;
    hl.num_wann= get<0>(wannier_data); //number of wannier functions
    const vector<string> hopping_lines= get<1>(wannier_data); //strings containing the hopping data
    for (auto line : hopping_lines){
        stringstream ss(line);
        ss.precision( numeric_limits<double>::digits10+2);
       
        hopping_list::cellID_t cellID; 
        for( auto& x : cellID)
            ss>>x;

        hopping_list::edge_t vertex_edge; 
        for( auto& x : vertex_edge){ ss>>x; x-=1; }//The input is assumed to be zero based        
        
        double re,im; ss>>re>>im; 
        hopping_list::value_t hop_value(re,im);

        if( sqrt( norm(hop_value) )> numeric_limits<double>::epsilon() )
            hl.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop_value,vertex_edge) } );
    }
return hl; 
};   

hopping_list wrap_in_supercell(const array<int, 3> cellDim,const hopping_list hl ){

    hopping_list sc_hl; //the hopping list for the supercell
    sc_hl.cellSizes=cellDim; 
    int numCells = 1; for(auto x: cellDim ) numCells*=x;
    sc_hl.num_wann = numCells*hl.num_wann;
    //Go through all the hopping list defined in the unit cell 
    //map the cell indexes into the new dimensions of the super cell
    //and add the values when it correspond
    for (auto const& key_hop : hl.hoppings){
        const auto hop = key_hop.second;
        auto cellID = get<0>(hop);
        auto value  = get<1>(hop);
        auto edge   = get<2>(hop);

        //wrap tag_indexes around the super_cell 
        for( size_t i=0; i < cellID.size(); i++)
            cellID[i]=( cellID[i]+cellDim[i])%cellDim[i];

        //compute the vertices shift
        int shift=1;
        for(auto x: cellID ) shift*=x*hl.num_wann;

        //Shift the edges appropiately    
        for( auto& vertice: edge )
            vertice+=shift;

        //look for the cell with the same ID and create or add the value
        string cellID_tag = get_tag(cellID,edge);
        if (sc_hl.hoppings.count(cellID_tag) ==  0)
            sc_hl.hoppings.insert( { cellID_tag, hopping_list::hopping_t(cellID,value,edge) } );
        else
            get<1>(sc_hl.hoppings[cellID_tag]) += value; 
    }
return sc_hl; 
};
