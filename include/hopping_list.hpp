#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>

using namespace std;

inline string tag_to_label(const array<int,3>& tag){
    string text_tag;
    for( auto& ti : tag)
        text_tag += to_string(ti)+" ";
return text_tag; 
}

struct hopping_list
{
    int num_wann;
    vector< tuple<int,int> > vertex_edge;
    map<string, tuple< array<int, 3>,complex<double> > > hoppings; 
};

hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  );

hopping_list wrap_in_supercell(const array<int, 3> cellDim, hopping_list hl );


inline int map_indexes3Dto1D( vector<int>& indexes, vector<int>& bounds  )
{
//    assert( indexes.size() == bounds.size() );
    int index = 0;
    for( int i = 0; i+1 < indexes.size(); i++ )
        index += ( (indexes[i]+bounds[i])%bounds[i] )*bounds[i+,1];
    index += ( indexes.back()+bounds.back() )%( bounds.back() )
    return index;
}

inline void save_hopping_list_as_csr(string output_filename, hopping_list& hl){


    for (auto const& hop : hl.hoppings){
        auto cellID_value_pair = hop.second;
        auto cellID = get<0>(cellID_value_pair);
        auto value = get<1>(cellID_value_pair);

        map_indexes3Dto1D
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
