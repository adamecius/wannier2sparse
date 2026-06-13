// Guards that save_supercell_as_csr (counting assembly, stable order) stays
// BYTE-IDENTICAL to a reference setFromTriplets assembly, for HAM and every
// operator family, including the maximal-duplicate collapse case {1,1,1}.
// This pins the deterministic accumulation order so the text CSR never drifts.
#include <vector>
#include <complex>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include "tbmodel.hpp"
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
using namespace std;

// Reference: the previous setFromTriplets-based assembly, kept here as the oracle.
static void ref_setFromTriplets(const hopping_list::cellID_t& cellDim,
                                const hopping_list& hl, const string& fname){
    const int WBB=hl.WannierBasisSize();
    const long Ntot=(long)cellDim[0]*cellDim[1]*cellDim[2];
    const size_t dim=(size_t)WBB*Ntot;
    SparseMatrix_t output(dim,dim);
    vector<Triplet_t> coeff; coeff.reserve(hl.hoppings.size()*Ntot);
    hopping_list::cellID_t cs;
    for(cs[2]=0;cs[2]<cellDim[2];++cs[2])for(cs[1]=0;cs[1]<cellDim[1];++cs[1])for(cs[0]=0;cs[0]<cellDim[0];++cs[0])
      for(size_t h=0;h<hl.hoppings.size();++h){
        hopping_list::cellID_t cid=get<0>(hl.hoppings[h]); auto v=get<1>(hl.hoppings[h]); auto e=get<2>(hl.hoppings[h]);
        for(int i=0;i<3;++i) cid[i]=((cid[i]+cs[i])+cellDim[i])%cellDim[i];
        int row=e[0]+index_aliasing(cs,cellDim)*WBB, col=e[1]+index_aliasing(cid,cellDim)*WBB;
        coeff.push_back(Triplet_t(row,col,v)); }
    output.setFromTriplets(coeff.begin(),coeff.end()); output.makeCompressed();
    ofstream f(fname.c_str());
    f<<dim<<" "<<output.nonZeros()<<endl;
    for(int k=0;k<output.outerSize();++k) for(SparseMatrix_t::InnerIterator it(output,k);it;++it) f<<it.value().real()<<" "<<it.value().imag()<<" "; f<<endl;
    for(int k=0;k<output.outerSize();++k) for(SparseMatrix_t::InnerIterator it(output,k);it;++it) f<<it.index()<<" "; f<<endl;
    for(int k=0;k<output.outerSize()+1;++k) f<<*(output.outerIndexPtr()+k)<<" "; f<<endl;
}
static string slurp(const string&p){ ifstream f(p); stringstream ss; ss<<f.rdbuf(); return ss.str(); }

int main(){
    tbmodel m; m.readOrbitalPositions("spin_graphene.xyz"); m.readUnitCell("spin_graphene.uc"); m.readWannierModel("spin_graphene_hr.dat");
    vector< pair<string,hopping_list> > ops;
    ops.push_back(make_pair(string("HAM"),  m.hl));
    ops.push_back(make_pair(string("VX"),   m.createHoppingCurrents_list(0)));
    ops.push_back(make_pair(string("SZ"),   m.createHoppingSpinDensity_list('z')));
    ops.push_back(make_pair(string("SY"),   m.createHoppingSpinDensity_list('y')));
    ops.push_back(make_pair(string("VYSZ"), m.createHoppingSpinCurrents_list(1,'z')));
    int n=0;
    hopping_list::cellID_t cds[4]={{1,1,1},{2,2,1},{3,3,1},{7,5,1}};
    for(int ci=0;ci<4;++ci) for(size_t oi=0;oi<ops.size();++oi){
        ref_setFromTriplets   (cds[ci], ops[oi].second, "ref.csr");
        save_supercell_as_csr (cds[ci], ops[oi].second, "got.csr");
        if(slurp("ref.csr")!=slurp("got.csr")){
            cerr<<"FAIL "<<ops[oi].first<<" at ("<<cds[ci][0]<<","<<cds[ci][1]<<","<<cds[ci][2]<<")\n"; return 1; }
        ++n;
    }
    cout<<"CSR COUNT-VS-SETFROMTRIPLETS EQUIVALENCE: PASSED ("<<n<<" cases, byte-identical)\n";
    return 0;
}
