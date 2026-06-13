#include "hopping_list.hpp"

hopping_list create_hopping_list( tuple<int, vector<int>, vector<string> > wannier_data  )
{
    hopping_list hl;
    const int num_wann = get<0>(wannier_data);
    hl.num_wann = num_wann;                                  //number of wannier functions
    const vector<int>&    ndegen        = get<1>(wannier_data); //WS degeneracy per R block
    const vector<string>& hopping_lines = get<2>(wannier_data); //strings with the hopping data

    // Hoppings arrive in blocks of num_wann^2 lines, one block per Wigner-Seitz
    // point, in the same order as ndegen. Each value is divided by its block's
    // degeneracy (the standard W90 ndegen normalization; a no-op when all are 1).
    // The block index is computed from the RAW line ordinal so the near-zero
    // filter below cannot desynchronize the counter.
    const long block_size = (long)num_wann * num_wann;
    long raw = 0;
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

        const long block = (block_size > 0) ? raw / block_size : 0;
        if( block < (long)ndegen.size() && ndegen[block] > 0 )
            hop_value /= (double)ndegen[block];
        ++raw;

        if( sqrt( norm(hop_value) )> numeric_limits<double>::epsilon() )
            hl.hoppings.push_back(hopping_list::hopping_t(cellID,hop_value,vertex_edge));
    }
return hl;
};

hopping_list wrap_in_supercell(const hopping_list::cellID_t& cellDim,const hopping_list hl ){

    hopping_list sc_hl;;  //the hopping list for the supercell
    sc_hl.SetWannierBasisSize(hl.WannierBasisSize()); //increase the supercell dimension;
    sc_hl.SetBounds(cellDim); //increase the supercell dimension;

    // Replicate the known Wannier connectivity into every translated cell.
    // The input already contains the neighbor graph, so the hot path should
    // only do integer wrapping and append generated hoppings.
    hopping_list::cellID_t cellShift;
    sc_hl.hoppings.reserve(hl.hoppings.size()*cellDim[0]*cellDim[1]*cellDim[2]);
    for(cellShift[2]=0; cellShift[2]< cellDim[2]; cellShift[2]++)
    for(cellShift[1]=0; cellShift[1]< cellDim[1]; cellShift[1]++)
    for(cellShift[0]=0; cellShift[0]< cellDim[0]; cellShift[0]++)
    {
        //Go through all the hopping list defined in the unit cell
        //map the cell indexes into the new dimensions of the super cell
        //and add the values when it correspond
        for (auto const& hop : hl.hoppings){
            auto cellID = get<0>(hop);
            auto value  = get<1>(hop);
            auto edge   = get<2>(hop);

             //Shift the cell given by the hopping vector
            //cellID = cellID + cellShift ;
            std::transform( cellID.begin(),cellID.end(),cellShift.begin(),cellID.begin(),  std::plus<int>() );

            //wrap tag_indexes around the super_cell
            for( size_t i=0; i < cellID.size(); i++)
                cellID[i]=( cellID[i]+cellDim[i])%cellDim[i];

            //Shift both the beggining and end of the edge.
            edge[0] += index_aliasing(cellShift,cellDim)*hl.WannierBasisSize(); //This is assume to be at the origin and shifted by cellShift
            edge[1] += index_aliasing(cellID,cellDim)*hl.WannierBasisSize(); //This is originally shidted by cellID and one added an aditional shift by cellShift.
            cellID = {0,0,0}; //After wrapping, all atoms belong to the supercell, therefore, cellID is always zero.

            // Append the hopping. If several unit-cell hoppings wrap onto the
            // same supercell edge, Eigen will combine them in setFromTriplets.
            assert( edge[0]< sc_hl.WannierBasisSize() );
            assert( edge[1]< sc_hl.WannierBasisSize() );
            assert( cellID== hopping_list::cellID_t({0,0,0}) );
            sc_hl.hoppings.push_back(hopping_list::hopping_t(cellID,value,edge));
        }
    }
    sc_hl.SetBounds(hopping_list::cellID_t({1,1,1})); //The wrapped cell is bounded in itself
return sc_hl;
};


void save_supercell_as_csr(const hopping_list::cellID_t& cellDim,
                           const hopping_list& hl, string output_filename)
{
    const int    WBB  = hl.WannierBasisSize();
    const long   Ntot = (long)cellDim[0]*cellDim[1]*cellDim[2];
    const size_t dim  = (size_t)WBB * Ntot;

    // Direct CSR (RowMajor) assembly by counting, bypassing setFromTriplets:
    //  - never holds a Triplet array and the Eigen matrix at the same time
    //    (~2.4x lower peak RSS, growing with size),
    //  - duplicate edges within a row are summed with a STABLE order that
    //    matches setFromTriplets, so the text CSR is byte-identical.
    // Output format is unchanged (see test/csr_count_equivalence).
    std::vector<long> rowptr(dim+1, 0);
    hopping_list::cellID_t cs;
    // pass 1: count entries per row (row depends only on cellShift and edge[0])
    for(cs[2]=0;cs[2]<cellDim[2];++cs[2])
    for(cs[1]=0;cs[1]<cellDim[1];++cs[1])
    for(cs[0]=0;cs[0]<cellDim[0];++cs[0])
        for(size_t h=0; h<hl.hoppings.size(); ++h){
            const auto edge = get<2>(hl.hoppings[h]);
            rowptr[ edge[0] + index_aliasing(cs,cellDim)*WBB + 1 ]++;
        }
    for(size_t r=0;r<dim;++r) rowptr[r+1]+=rowptr[r];

    const long nnz_raw = rowptr[dim];
    std::vector<int>                  col(nnz_raw);
    std::vector<std::complex<double>> val(nnz_raw);
    std::vector<long>                 pos(rowptr.begin(), rowptr.end()-1);
    // pass 2: scatter into per-row segments, in the same loop order as before
    for(cs[2]=0;cs[2]<cellDim[2];++cs[2])
    for(cs[1]=0;cs[1]<cellDim[1];++cs[1])
    for(cs[0]=0;cs[0]<cellDim[0];++cs[0])
        for(size_t h=0; h<hl.hoppings.size(); ++h){
            hopping_list::cellID_t cellID = get<0>(hl.hoppings[h]);
            const auto value = get<1>(hl.hoppings[h]);
            const auto edge  = get<2>(hl.hoppings[h]);
            for(int i=0;i<3;++i) cellID[i]=((cellID[i]+cs[i])+cellDim[i])%cellDim[i];
            const int row = edge[0] + index_aliasing(cs,    cellDim)*WBB;
            const int c   = edge[1] + index_aliasing(cellID,cellDim)*WBB;
            assert(row<(int)dim); assert(c<(int)dim);
            const long p = pos[row]++; col[p]=c; val[p]=value;
        }

    // compress per row IN PLACE: stable sort by column + sum duplicates.
    // write cursor never overtakes the read cursor (merged <= raw per row).
    std::vector<long> rp(dim+1, 0);
    std::vector< std::pair<int,std::complex<double> > > tmp;
    for(size_t r=0;r<dim;++r){
        const long s=rowptr[r], e=rowptr[r+1];
        tmp.clear(); tmp.reserve(e-s);
        for(long k=s;k<e;++k) tmp.push_back( std::make_pair(col[k],val[k]) );
        std::stable_sort(tmp.begin(),tmp.end(),
            [](const std::pair<int,std::complex<double> >&a,
               const std::pair<int,std::complex<double> >&b){return a.first<b.first;});
        long w = rp[r]; size_t k=0;
        while(k<tmp.size()){
            int cc=tmp[k].first; std::complex<double> acc=tmp[k].second; size_t j=k+1;
            while(j<tmp.size() && tmp[j].first==cc){ acc+=tmp[j].second; ++j; }
            col[w]=cc; val[w]=acc; ++w; k=j;
        }
        rp[r+1]=w;
    }
    const long nnz = rp[dim];

    std::ofstream matrix_file ( output_filename.c_str() );
    matrix_file<<dim<<" "<<nnz<<std::endl;
    for(long k=0;k<nnz;++k) matrix_file<<val[k].real()<<" "<<val[k].imag()<<" ";
    matrix_file<<std::endl;
    for(long k=0;k<nnz;++k) matrix_file<<col[k]<<" ";
    matrix_file<<std::endl;
    for(size_t r=0;r<=dim;++r) matrix_file<<rp[r]<<" ";
    matrix_file<<std::endl;
    matrix_file.close();
return ;
}
