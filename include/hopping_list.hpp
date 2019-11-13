#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>
#include <iostream>
#include <limits>

using namespace std;

struct hopping_list
{
    typedef complex<double> value_t;
    typedef array<int, 2> edge_t;
    typedef array<int, 3> cellID_t;
    typedef tuple< cellID_t,value_t,edge_t > hopping_t;

    hopping_list():cellSizes({1,1,1}), num_wann(0){};

    cellID_t Bounds()
    {
        return array<int, 3>({cellSizes[0],cellSizes[1],cellSizes[2]});
    };

    int cellID_index(const cellID_t cidx )
    {
        const cellID_t bounds = this->Bounds();
       
        int index = 0;
        for( int i = 0; i+1 < cidx.size(); i++ )
            index += ( (cidx[i]+bounds[i])%bounds[i] )*bounds[i+1];
        index += ( cidx.back()+bounds.back() )%( bounds.back() );
        return index;
    };
    bool operator ==(hopping_list& y )
    {
        bool list_equal = true;
        for( auto const& elem: y.hoppings )
        {
            auto key = elem.first;
            if( this->hoppings.count(key) == 0 )
            {
                list_equal = false;    
                break;
            }
            list_equal*= (bool)(get<0>(this->hoppings[key])==get<0>(elem.second ));
//          std::cou<<"FIX_EQUAL"  list_equal*=  (bool)( sqrt(norm( get<1>(this->hoppings[key]) - get<1>(elem.second) )) < numeric_limits<double>::epsilon() );
            list_equal*= (bool)(get<2>(this->hoppings[key])==get<2>(elem.second ));
        }
        return  ( this->num_wann ==y.num_wann)&&
                (this->cellSizes==y.cellSizes)&&
                list_equal;
    }

    int num_wann;
    cellID_t cellSizes;
    map<string, hopping_t > hoppings; 
};

hopping_list create_hopping_list( tuple<int, vector<string> > wannier_data  );

hopping_list wrap_in_supercell(const array<int, 3> cellDim, hopping_list hl );

inline string get_tag(const hopping_list::cellID_t& cid,const hopping_list::edge_t edge){
    string text_tag;
    for( auto& ti : cid)
        text_tag += to_string(ti)+" ";
    for( auto& ti : edge)
        text_tag += to_string(ti)+" ";
return text_tag; 
}

inline array<int,5> tag_to_indices(const string& tag){
    array<int,5> indices;
    stringstream ss(tag);
    for( auto& ti : indices)
        ss>>ti;
return indices; 
}
#include <iostream>



inline void save_hopping_list_as_csr(string output_filename, hopping_list& hl){

    hopping_list sc_hl;
    
    for(auto const& key_hop : hl.hoppings){

        const auto key = key_hop.first;
        const auto hop = key_hop.second;
        auto cellID = get<0>(hop);
        auto value  = get<1>(hop);
        auto edge   = get<2>(hop);
        
        int 
        row=hl.cellID_index(cellID)*hl.num_wann + get<0>(edge),
        col=hl.cellID_index(cellID)*hl.num_wann + get<1>(edge);

        std::cout<<row<<" "<<col<<" "<<value<<std::endl;

    }
return ; 
}