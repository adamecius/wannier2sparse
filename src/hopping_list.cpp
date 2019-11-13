#include "hopping_list.hpp"
#include <functional>
#include<iostream>
#include<limits>



hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  )
{
    typedef complex<double> complexd;
    typedef tuple< array<int, 3>,complexd> hopping_value;
    hopping_list hl;

    hl.num_wann= get<0>(wannier_data); //number of wannier functions
    const vector<string> hopping_lines= get<1>(wannier_data); //strings containing the hopping data

    for (auto line : hopping_lines){
        stringstream ss(line);
               
        array<int, 3> tag; string text_tag; 
        for( auto& ti : tag)
            ss>>ti;
        text_tag = tag_to_label(tag);

        int forb,torb; 
        ss>>forb>>torb; //initial and final vertex defining a hopping

        double rv,im; 
        ss>>rv>>im;
        complexd hop_value(rv,im);

        if( sqrt( std::norm(hop_value) )> numeric_limits<double>::epsilon() )
        {
            hl.hoppings.insert( {text_tag, hopping_value(tag,complexd(rv,im) ) } );
            hl.vertex_edge.push_back(tuple<int,int>(forb,torb) );
        }
    }

return hl; 
}   

hopping_list wrap_in_supercell(const array<int, 3> cellDim, hopping_list hl ){

    typedef complex<double> complexd;
    typedef tuple< array<int, 3>,complexd> hopping_value;
    hopping_list sc_hl; //the hopping list for the supercell

    //Go through all the hopping list defined in the unit cell 
    //map the cell indexes into the new dimensions of the super cell
    //and add the values when it correspond
    for (auto const& hop : hl.hoppings){
        auto cellID_value_pair = hop.second;
        auto cellID = get<0>(cellID_value_pair);
        auto value = get<1>(cellID_value_pair);

        //wrap tag_indexes around the super_cell 
        for( size_t i=0; i < cellID.size(); i++)
            cellID[i]=( cellID[i]+cellDim[i])%cellDim[i];
        
        //look for the cell with the same ID and create or add the value
        string cellID_tag = tag_to_label(cellID);
        if (sc_hl.hoppings.count(cellID_tag) ==  0)
            sc_hl.hoppings.insert( { cellID_tag, hopping_value(cellID,value) } );
        else
            get<1>(sc_hl.hoppings[cellID_tag]) += value; 
    }
return sc_hl; 
}

