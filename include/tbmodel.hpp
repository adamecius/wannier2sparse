#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>
#include <iostream>
#include <limits>
#include <cassert>
#include <functional>
#include<iostream>
#include<limits>
#include<algorithm>
#include "wannier_parser.hpp"

using namespace std;


double volume( const array < array<double,3> , 3 >&  uc )
{
    return  std::fabs(
            ( uc[0][1]*uc[1][2]-uc[1][1]*uc[0][2])*uc[2][0]+
            ( uc[0][2]*uc[1][0]-uc[0][0]*uc[1][2])*uc[2][1]+
            ( uc[0][0]*uc[1][1]-uc[0][1]*uc[1][0])*uc[2][2]);
}

class tbmodel
{
    public:
    typedef tuple<string, array<double, 3> > orbPos_t;
    typedef array < array<double,3> , 3 > unitCell_t;


    inline void readUnitCell(const string inputfile)
    {
        lat_vec= read_unit_cell_file("graphene.uc");

    }

    inline void readOrbitalPositions(const string inputfile)
    {
        orbPos_list= read_xyz_file("graphene.xyz");
    };

    inline void readWannierModel(const string inputfile)
    {
        hl = create_hopping_list(read_wannier_file(inputfile));
    };

    hopping_list createHoppingCurrents_list(const int dir)
    {
        assert( dir <3 && dir >=0 ); 
        assert(volume(lat_vec) > 0 );
        assert(orbPos_list.size()==hl.WannierBasisSize());

        hopping_list chl = this->hl ;
        for( auto& elem: chl.hoppings )
        {
            auto  key = elem.first;
            auto  tag  = get<0>(elem.second);
            auto  edge = get<2>(elem.second);
            auto  val  = get<1>(elem.second);
            array<double,3>
            pos( 
                {
                    ( get<1>(orbPos_list[edge[0]])[0] - get<1>(orbPos_list[edge[1]])[0] ) + ( tag[0]*lat_vec[0][0] + tag[1]*lat_vec[1][0] + tag[2]*lat_vec[2][0] ),
                    ( get<1>(orbPos_list[edge[0]])[1] - get<1>(orbPos_list[edge[1]])[1] ) + ( tag[0]*lat_vec[0][1] + tag[1]*lat_vec[1][1] + tag[2]*lat_vec[2][1] ),
                    ( get<1>(orbPos_list[edge[0]])[2] - get<1>(orbPos_list[edge[1]])[2] ) + ( tag[0]*lat_vec[0][2] + tag[1]*lat_vec[1][2] + tag[2]*lat_vec[2][2] )  
                }                     
                );
            get<1>(elem.second) *=  hopping_list::value_t(0.0,pos[dir]);                
            std::cout<<key<<" | "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
        }
        return chl;
    };


    hopping_list createHoppingDensity_list(const int dir)
    {
        assert( dir <3 && dir >=0 ); 
        assert(volume(lat_vec) > 0 );
        assert(orbPos_list.size()==hl.WannierBasisSize());

        hopping_list chl = this->hl ;
        for( auto& elem: chl.hoppings )
        {
            auto  key = elem.first;
            auto  tag  = get<0>(elem.second);
            auto  edge = get<2>(elem.second);
            auto  val  = get<1>(elem.second);
            if( tag[0]== 0 && tag[1]== 0 && tag[2]== 0 )
            {
                array<double,3>
                pos( 
                    {
                        ( get<1>(orbPos_list[edge[0]])[0] - get<1>(orbPos_list[edge[1]])[0] ) ,
                        ( get<1>(orbPos_list[edge[0]])[1] - get<1>(orbPos_list[edge[1]])[1] ) ,
                        ( get<1>(orbPos_list[edge[0]])[2] - get<1>(orbPos_list[edge[1]])[2] )   
                    }                     
                    );
                get<1>(elem.second) *=  hopping_list::value_t(0.0,pos[dir]);                
                std::cout<<key<<" | "<<pos[0]<<" "<<pos[1]<<" "<<pos[2]<<std::endl;
            }
        }
        return chl;
    };


    int num_orbs;
    unitCell_t lat_vec;
    vector< orbPos_t > orbPos_list;
    hopping_list hl;
};
