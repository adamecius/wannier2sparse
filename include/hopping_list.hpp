/**
 * @file hopping_list.hpp
 * @brief Hopping list: the core data structure for Wannier90 real-space operators.
 *
 * A hopping_list stores the primitive-cell operator \f$O_{ij}(R)\f$ as a flat list
 * of (cellID, value, edge) tuples. The storage is append-only and intentionally
 * preserves duplicate edges; they are summed later during sparse-matrix assembly.
 */
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
#include"wannier_parser.hpp"   // wsvec_entry (lightweight; no Eigen)
#include <iostream>

using namespace std;

/**
 * @brief Flatten a 3D cell index to a linear cell index under periodic wrapping.
 * @param index cell index (may be negative)
 * @param bound supercell dimensions (must be > 0)
 * @return linear index in [0, bound[0]*bound[1]*bound[2])
 */
inline int
index_aliasing(const array<int, 3>& index,const array<int, 3>& bound )
{
    return ( (index[2]+bound[2])%bound[2] * bound[1] + (index[1]+bound[1])%bound[1] ) * bound[0] + (index[0]+bound[0])%bound[0] ;
}


/**
 * @brief Container for a Wannier90 real-space operator.
 *
 * Stores \f$O_{ij}(R)\f$ as a flat append-only list. The old map keyed by a string
 * made from (cellID, edge) was replaced by this vector-backed storage because the
 * key was never used for lookup in the production pipeline and forced a string
 * allocation plus tree insertion for every replicated hopping.
 */
struct hopping_list
{
    typedef complex<double> value_t;            ///< Complex hopping value (eV, eV·Å, ħ/2, ...)
    typedef array<int, 2> edge_t;               ///< (initial orbital, final orbital), zero-based
    typedef array<int, 3> cellID_t;             ///< (Rx, Ry, Rz) Wigner-Seitz cell vector
    typedef tuple< cellID_t,value_t,edge_t > hopping_t; ///< Single hopping record

    /**
     * @brief Storage adapter that keeps a vector<hopping_t> while preserving a
     *        legacy string-keyed insert/lookup interface for tests.
     */
    struct hopping_storage : public vector<hopping_t>
    {
        using vector<hopping_t>::insert;
        using vector<hopping_t>::operator[];

        /**
         * @brief Append a tagged hopping, ignoring the tag.
         * @param tagged_hop pair whose second element is appended
         */
        void insert(const pair<string, hopping_t>& tagged_hop)
        {
            this->push_back(tagged_hop.second);
        }

        /**
         * @brief Linear search for a hopping matching the (cellID, edge) tag.
         * @param tag string tag produced by make_tag()
         * @return reference to the first matching hopping
         * @warning Aborts with assert if no match is found.
         */
        hopping_t& operator[](const string& tag)
        {
            for( auto& hop : *this )
                if( make_tag(hop) == tag )
                    return hop;

            std::cout<<"The key: "<<tag<<" was not found in hopping_storage"<<std::endl;
            assert(false);
            return this->front();
        }

        /**
         * @brief Const linear search for a hopping matching the tag.
         * @param tag string tag produced by make_tag()
         * @return const reference to the first matching hopping
         * @warning Aborts with assert if no match is found.
         */
        const hopping_t& operator[](const string& tag) const
        {
            for( auto const& hop : *this )
                if( make_tag(hop) == tag )
                    return hop;

            std::cout<<"The key: "<<tag<<" was not found in hopping_storage"<<std::endl;
            assert(false);
            return this->front();
        }

        /**
         * @brief Build the legacy textual tag for a hopping.
         * @param hop hopping tuple
         * @return "Rx Ry Rz i j " string
         */
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

    /**
     * @brief Default constructor: one unit cell, zero Wannier functions.
     */
    hopping_list():cellSizes({1,1,1}), num_wann(0){};

    /**
     * @brief Number of Wannier orbitals in the current bounding box.
     * @return num_wann
     */
    inline int WannierBasisSize() const
    {
        return this-> num_wann ;
    };

    /**
     * @brief Set the number of Wannier orbitals per unit cell.
     * @param num_wann must be > 0
     */
    inline void SetWannierBasisSize(const int num_wann)
    {
        assert( num_wann > 0 );
        this-> num_wann  = num_wann;
        return ;
    };

    /**
     * @brief Set the bounding box (number of unit cells) and scale num_wann.
     * @param cellSizes dimensions along each lattice vector (each > 0)
     *
     * num_wann is multiplied by the product of the three dimensions, so after
     * SetBounds the object describes the full supercell orbital count.
     */
    inline void SetBounds(const cellID_t& cellSizes)
    {
        assert( this-> num_wann > 0 && cellSizes[0]>0&& cellSizes[1]>0&&cellSizes[2]>0 );
        this->cellSizes = cellSizes;
        for(const auto& x: cellSizes )
            this->num_wann *=x;
        return ;
    };

    /**
     * @brief Current bounding box.
     * @return (Nx, Ny, Nz)
     */
    cellID_t Bounds()
    {
        return array<int, 3>({cellSizes[0],cellSizes[1],cellSizes[2]});
    };

    /**
     * @brief Linear cell index within the current bounding box.
     * @param cidx cell index (may be negative; wrapped)
     * @return linear index
     */
    int cellID_index(const cellID_t cidx )
    {
        const cellID_t bounds = this->Bounds();

        int index = 0;
        for( int i = 0; i+1 < cidx.size(); i++ )
            index += ( (cidx[i]+bounds[i])%bounds[i] )*bounds[i+1];
        index += ( cidx.back()+bounds.back() )%( bounds.back() );
        return index;
    };

    /**
     * @brief Compare two hopping lists for approximate equality.
     * @param y other hopping list
     * @return true if size, bounds, and all summed hoppings agree
     *
     * Hoppings are grouped by (cellID, edge), summed, and compared with a
     * relative tolerance of machine epsilon.
     */
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

    int num_wann;        ///< number of Wannier orbitals in the current bounding box
    cellID_t cellSizes;  ///< current bounding box dimensions

    /**
     * @brief Flat append-only hopping storage.
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
 * @brief Parse Wannier90 `_hr.dat` data into a hopping list.
 * @param wannier_data tuple of (num_wann, ndegen, hopping_lines)
 * @return populated hopping_list
 *
 * Each hopping is divided by the Wigner-Seitz degeneracy of its R-block (standard
 * W90 normalization; a no-op when every ndegen is 1). Orbital indices are
 * converted from Wannier's one-based convention to zero-based indices used
 * internally.
 */
hopping_list create_hopping_list( tuple<int, vector<int>, vector<string> > wannier_data  );

/**
 * @brief Apply the Wigner-Seitz minimum-image correction.
 * @param hl input hopping list
 * @param wsvec correction data from `_wsvec.dat`
 * @return corrected hopping list
 *
 * Each hopping (R, i, j) is replaced by its T.size() images at R + T, each
 * weighted 1/T.size(); this composes with the ndegen division already done in
 * create_hopping_list. Hoppings with no matching wsvec record, and an empty
 * wsvec, are passed through unchanged (so it is a no-op without use_ws_distance).
 */
hopping_list apply_wsvec(const hopping_list& hl, const vector<wsvec_entry>& wsvec);

/**
 * @brief Replicate a unit-cell hopping list over a periodic supercell.
 * @param cellDim supercell dimensions (each >= 1)
 * @param hl input hopping list
 * @return expanded hopping list
 *
 * Each original hopping is translated to every cell in cellDim, wrapped back
 * into the supercell, and converted to supercell orbital indices. The resulting
 * list may contain duplicate edges; those are cheaper to append here and are
 * later combined during sparse matrix assembly.
 */
hopping_list wrap_in_supercell(const hopping_list::cellID_t& cellDim,const hopping_list hl);

/**
 * @brief Fused single-pass supercell expansion + CSR export.
 * @param cellDim supercell dimensions
 * @param hl primitive-cell hopping list
 * @param output_filename path for the `.CSR` output
 *
 * Replicates and PBC-wraps the primitive-cell hoppings straight into flat
 * (row,col,value) triplets and writes the CSR file, without ever materialising
 * the intermediate supercell hopping_list. The output is byte-identical to the
 * two-stage save_hopping_list_as_csr(output, wrap_in_supercell(cellDim, hl)),
 * while avoiding the extra full-supercell array. This is the production export
 * path.
 */
void save_supercell_as_csr(const hopping_list::cellID_t& cellDim,
                           const hopping_list& hl, string output_filename);

/**
 * @brief Detect minimum-image aliasing collisions.
 * @param hl primitive-cell hopping list
 * @param cellDim supercell dimensions
 * @return human-readable messages, one per colliding pair
 *
 * Two distinct R for the same (i,j) that collapse onto the same supercell bond
 * R mod cellDim are reported. Empty result means the supercell is safe.
 */
inline std::vector<std::string>
aliasing_collisions(const hopping_list& hl, const hopping_list::cellID_t& cellDim)
{
    auto wrap = [&](int v, int n){ return ((v % n) + n) % n; };
    std::map<std::array<int,5>, std::array<int,3> > seen;
    std::vector<std::string> out;
    for (const auto& h : hl.hoppings)
    {
        const auto R = std::get<0>(h); const auto e = std::get<2>(h);
        const std::array<int,5> wk = { wrap(R[0],cellDim[0]), wrap(R[1],cellDim[1]),
                                       wrap(R[2],cellDim[2]), e[0], e[1] };
        auto it = seen.find(wk);
        if (it == seen.end()) seen[wk] = {R[0],R[1],R[2]};
        else if (it->second != std::array<int,3>{R[0],R[1],R[2]})
            out.push_back("(i,j)=(" + std::to_string(e[0]) + "," + std::to_string(e[1]) +
                          "): R=(" + std::to_string(it->second[0]) + "," + std::to_string(it->second[1]) +
                          "," + std::to_string(it->second[2]) + ") and R=(" + std::to_string(R[0]) +
                          "," + std::to_string(R[1]) + "," + std::to_string(R[2]) + ")");
    }
    return out;
}

/**
 * @brief Abort if expanding hl into cellDim would alias distinct R onto the same bond.
 * @param hl primitive-cell hopping list
 * @param cellDim supercell dimensions
 *
 * Real-space minimum image needs N >= 2*range+1 per axis. Used by the CLI before
 * export. wrap_in_supercell/save_supercell_as_csr themselves do NOT call this;
 * they permit folding (the equivalence/folding unit tests rely on small-N
 * collapse).
 */
void guard_minimum_image(const hopping_list& hl, const hopping_list::cellID_t& cellDim);

/**
 * @brief Assemble the expanded supercell of hl as an Eigen sparse matrix.
 * @param cellDim supercell dimensions
 * @param hl primitive-cell hopping list
 * @return sparse matrix of the expanded supercell operator
 *
 * Replicates, PBC-wraps, and sums duplicate edges. Lets operators be combined
 * algebraically after expansion (e.g. the spin current J = 1/2{V,S}).
 */
SparseMatrix_t supercell_matrix(const hopping_list::cellID_t& cellDim,
                                const hopping_list& hl);

/**
 * @brief Write a sparse matrix to the CSR text format.
 * @param matrix sparse matrix to write
 * @param output_filename output path
 *
 * Format: dim nnz / real imag ... / column indices ... / row pointers ...
 * Factored out so both hopping-list export and derived (matrix) operators share
 * exactly one writer.
 */
void write_csr(const SparseMatrix_t& matrix, string output_filename);

/**
 * @brief Build the legacy textual tag for a hopping.
 * @param cid cellID
 * @param edge orbital edge
 * @return "Rx Ry Rz i j " string
 *
 * Retained for tests and diagnostics. Not used as the primary storage key.
 */
inline string get_tag(const hopping_list::cellID_t& cid,const hopping_list::edge_t edge){
    string text_tag;
    for( auto& ti : cid)
        text_tag += to_string(ti)+" ";
    for( auto& ti : edge)
        text_tag += to_string(ti)+" ";
return text_tag;
}

/**
 * @brief Parse a legacy textual tag back into indices.
 * @param tag tag string
 * @return (Rx, Ry, Rz, i, j)
 */
inline array<int,5> tag_to_indices(const string& tag){
    array<int,5> indices;
    stringstream ss(tag);
    for( auto& ti : indices)
        ss>>ti;
return indices;
}

/**
 * @brief Legacy two-stage export: expand to a hopping_list, then write CSR.
 * @param output_filename output path
 * @param hl hopping list to export
 *
 * Prefer save_supercell_as_csr for large supercells; this function is kept for
 * tests and small cases.
 */
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
