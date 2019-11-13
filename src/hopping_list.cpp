#include "hopping_list.hpp"

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
        for( auto& ti : tag){
            ss>>ti;
            text_tag += to_string(ti)+" ";
        }

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

const int hl.num_wann=10;
const array<int, 3> sc_dims( {2,2,2} ) ; 
const int num_cells = sc_dims[0]*sc_dims[1]*sc_dims[2];
const vector< array<int, 3> > sc_tags(num_cells ) ; 
vector< tuple<int,int> > uc_vertex_edge( hl.vertex_edge);

for (auto const& sc_tag: sc_tags ){
    vertexExtender expanded_vertices(sc_tag,uc_vertex_edge);
    hl.vertex_edge.append( expanded_vertices.begin(),expanded_vertices.end() )
}

return hl; }
/*
for (auto const& hop : hl.hoppings){
    
     array<int, 3> tag = get<0>(hop.second);
    cout << hop.first<< endl;
    cout << hop.first  // string (key)
         << ':'
         << tag[0]<<" "<<tag[1]<<" "<<tag[2]<<endl
         << get<1>(hop.second) // string's value
         << endl ;
}    */
