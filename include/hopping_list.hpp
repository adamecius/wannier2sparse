#ifndef HOPPING_LIST
#define HOPPING_LIST
#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>
#include <fstream>
#include <iostream>
#include <limits>
#include <cassert>
#include <functional>
#include<iostream>
#include<limits>
#include<algorithm>
#include"sparse_matrix.hpp"
#include <iostream>

using namespace std;

inline int
index_aliasing(const array<int, 3>& index,const array<int, 3>& bound )
{
    return ( (index[2]+bound[2])%bound[2] * bound[1] + (index[1]+bound[1])%bound[1] ) * bound[0] + (index[0]+bound[0])%bound[0] ;
}


struct hopping_list
{
    typedef complex<double> value_t;
    typedef array<int, 2> edge_t;
    typedef array<int, 3> cellID_t;
    typedef tuple< cellID_t,value_t,edge_t > hopping_t;

    struct hopping_storage : public vector<hopping_t>
    {
        using vector<hopping_t>::insert;
        using vector<hopping_t>::operator[];

        void insert(const pair<string, hopping_t>& tagged_hop)
        {
            this->push_back(tagged_hop.second);
        }

        hopping_t& operator[](const string& tag)
        {
            for( auto& hop : *this )
                if( make_tag(hop) == tag )
                    return hop;

            std::cout<<"The key: "<<tag<<" was not found in hopping_storage"<<std::endl;
            assert(false);
            return this->front();
        }

        const hopping_t& operator[](const string& tag) const
        {
            for( auto const& hop : *this )
                if( make_tag(hop) == tag )
                    return hop;

            std::cout<<"The key: "<<tag<<" was not found in hopping_storage"<<std::endl;
            assert(false);
            return this->front();
        }

        static string make_tag(const hopping_t& hop)
        {
            string tag;
            for( auto x : get<0>(hop) )
                tag += to_string(x)+" ";
            for( auto x : get<2>(hop) )
                tag += to_string(x)+" ";
            return tag;
        }
    };

    hopping_list():cellSizes({1,1,1}), num_wann(0){};

    inline int WannierBasisSize() const
    {
        return this-> num_wann ;
    };

    inline void SetWannierBasisSize(const int num_wann)
    {
        assert( num_wann > 0 );
        this-> num_wann  = num_wann;
        return ;
    };

    inline void SetBounds(const cellID_t& cellSizes)
    {
        assert( this-> num_wann > 0 && cellSizes[0]>0&& cellSizes[1]>0&&cellSizes[2]>0 );
        this->cellSizes = cellSizes;
        for(const auto& x: cellSizes )
            this->num_wann *=x;
        return ;
    };

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
        auto reduce_by_tag = [](const hopping_storage& hoppings)
        {
            map<string, hopping_t> reduced;
            for( auto const& hop : hoppings )
            {
                const auto tag = hopping_storage::make_tag(hop);
                auto match = reduced.find(tag);
                if( match == reduced.end() )
                    reduced.insert({tag, hop});
                else
                    get<1>(match->second) += get<1>(hop);
            }
            return reduced;
        };

        const auto x_hoppings = reduce_by_tag(this->hoppings);
        const auto y_hoppings = reduce_by_tag(y.hoppings);
        bool list_equal = true;
        for( auto const& elem: y_hoppings )
        {
            auto key = elem.first;
            if( x_hoppings.count(key) == 0 )
            {
                std::cout<<"The key: "<<key<<" was not found when comparing the hopping lists"<<std::endl;
                list_equal = false;
                break;
            }
            const auto& x_hop = x_hoppings.at(key);
            list_equal*= (bool)(get<0>(x_hop)==get<0>(elem.second));
            list_equal*= (bool)(get<2>(x_hop)==get<2>(elem.second));

            auto val_diff= (get<1>(x_hop)-get<1>(elem.second ))/2.0;
            auto val_sum = (get<1>(x_hop)+get<1>(elem.second ))/2.0;

            if(!( val_diff.real()==0&& val_diff.imag()==0 ) )
            {
                list_equal*= (bool)(
                                std::fabs(val_diff.real()/val_sum.real()) < std::numeric_limits<double>::epsilon() &&
                                std::fabs(val_diff.imag()/val_sum.imag()) < std::numeric_limits<double>::epsilon()
                                );
                if(!list_equal)
                    std::cout<<"The keys "<<key<<" have hoppings with a percentile difference higher than "<< std::numeric_limits<double>::epsilon()<<std::endl;
            }
        }
        return  ( this->num_wann ==y.num_wann)&&
                (this->cellSizes==y.cellSizes)&&
                (x_hoppings.size()==y_hoppings.size())&&
                list_equal;
    }

    int num_wann;
    cellID_t cellSizes;
    /**
     * Flat append-only hopping storage.
     *
     * Older versions used a map keyed by a string made from `(cellID, edge)`.
     * The key was never used for lookup in the production pipeline; it only
     * forced a string allocation and tree insertion for every replicated
     * hopping. Keeping the terms in insertion order makes supercell wrapping
     * linear in the number of generated hoppings. Duplicate `(row, col)` terms
     * are intentionally preserved here and are summed by Eigen when CSR
     * triplets are assembled.
     */
    hopping_storage hoppings;
};

/**
 * Parse Wannier90 `_hr.dat` data into a hopping list.
 *
 * Takes (num_wann, ndegen, hopping_lines). Each hopping is divided by the
 * Wigner-Seitz degeneracy of its R-block (standard W90 normalization; a no-op
 * when every ndegen is 1). The parser keeps all nonzero terms from the
 * unit-cell Hamiltonian. Orbital indices are converted from Wannier's one-based
 * convention to zero-based indices used internally.
 */
hopping_list create_hopping_list( tuple<int, vector<int>, vector<string> > wannier_data  );

/**
 * Replicate a unit-cell hopping list over a periodic supercell.
 *
 * Each original hopping is translated to every cell in `cellDim`, wrapped back
 * into the supercell, and converted to supercell orbital indices. The resulting
 * list may contain duplicate edges; those are cheaper to append here and are
 * later combined during sparse matrix assembly.
 */
hopping_list wrap_in_supercell(const hopping_list::cellID_t& cellDim,const hopping_list hl);

/**
 * Fused single-pass supercell expansion + CSR export.
 *
 * Replicates and PBC-wraps the primitive-cell hoppings straight into flat
 * (row,col,value) triplets and writes the CSR file, without ever materialising
 * the intermediate supercell hopping_list.  Triplets are emitted in the same
 * loop order as wrap_in_supercell, so setFromTriplets combines duplicate edges
 * in the same order: the output is byte-identical to the two-stage
 * save_hopping_list_as_csr(output, wrap_in_supercell(cellDim, hl)), while
 * avoiding the extra full-supercell array.  This is the production export path.
 */
void save_supercell_as_csr(const hopping_list::cellID_t& cellDim,
                           const hopping_list& hl, string output_filename);

// Assemble the expanded supercell of `hl` as an Eigen sparse matrix (replicate +
// PBC-wrap, then sum duplicate edges). Lets operators be combined algebraically
// after expansion (e.g. the spin current J = 1/2{V,S}).
SparseMatrix_t supercell_matrix(const hopping_list::cellID_t& cellDim,
                                const hopping_list& hl);

// Write a sparse matrix to the CSR text format used throughout the tool:
//   dim nnz / real imag ... / column indices ... / row pointers ...
// Factored out so both hopping-list export and derived (matrix) operators share
// exactly one writer.
void write_csr(const SparseMatrix_t& matrix, string output_filename);

/**
 * Build the legacy textual tag for a hopping.
 *
 * This is retained for tests and diagnostics. It is not used as the primary
 * storage key anymore because tag construction was the hot-path bottleneck in
 * large supercell generation.
 */
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

inline void save_hopping_list_as_csr(string output_filename,const hopping_list& hl)
{
    const size_t dim = hl.WannierBasisSize();
    SparseMatrix_t output(dim,dim);
    std::vector<Triplet_t> coefficients;            // list of non-zeros coefficients
    coefficients.reserve(hl.hoppings.size());
    for(auto const& elem : hl.hoppings)
    {
        const auto value  = get<1>(elem);
        const auto edge   = get<2>(elem);
        coefficients.push_back(Triplet_t(edge[0],edge[1],value) );
    }
    output.setFromTriplets(coefficients.begin(), coefficients.end());
    output.makeCompressed();


	std::ofstream matrix_file ( output_filename.c_str()) ;

	//READ DIMENSION OF THE MATRIX
	matrix_file<<dim<<" "<<output.nonZeros()<<std::endl;

    //save values first
    for (int k=0; k<output.outerSize(); ++k)
    for (SparseMatrix_t::InnerIterator it(output,k); it; ++it)
        matrix_file<<it.value().real()<<" "<<it.value().imag()<<" ";
    matrix_file<<std::endl;

    //save the columns
    for (int k=0; k<output.outerSize(); ++k)
    for (SparseMatrix_t::InnerIterator it(output,k); it; ++it)
        matrix_file<<it.index()<<" ";
    matrix_file<<std::endl;

    //save the indices to columns
    for (int k=0; k<output.outerSize()+1; ++k)
        matrix_file<<*( output.outerIndexPtr() + k ) <<" ";
    matrix_file<<std::endl;

    matrix_file.close();

return ;
}

#endif
