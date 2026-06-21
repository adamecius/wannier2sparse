/**
 * @file tbmodel.hpp
 * @brief Tight-binding model and operator builders.
 *
 * Reads the Wannier90 inputs (lattice vectors, orbital positions, Hamiltonian) and
 * builds derived real-space operators: velocity/current, density, spin density, and
 * spin current. Spin operators rely on orbital labels marking the spin channel
 * (`_s+_` / `_s-_`); spinless models produce zero for those operators.
 */
#ifndef TBMODELS
#define TBMODELS
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
#include <stdexcept>
#include <functional>
#include <cmath>
#include<iostream>
#include<limits>
#include<algorithm>
#include "wannier_parser.hpp"
#include "hopping_list.hpp"

using namespace std;


/**
 * @brief Print and assert equality.
 * @tparam T comparable type
 * @param x first value
 * @param y second value
 */
template<typename T>
void assert_equal(const T x,const T y )
{
    if( x != y )
        std::cerr<<"ASSERT_EQUAL FAILED: "<<x<<" != "<<y<<std::endl;
    assert(x==y);
    return ;
}

/**
 * @brief Convert a cell vector to Cartesian coordinates.
 * @param tag integer cell index (Rx,Ry,Rz)
 * @param lat_vecs lattice vectors (rows)
 * @return Cartesian vector (Angstrom)
 */
array<double,3> tag2cartesian(const array<int,3>& tag,const array < array<double,3> , 3 >& lat_vecs  )
{
    array<double,3> cart_vect = {0,0,0};
    for(int i =0 ; i < cart_vect.size(); i++ )
    for(int j =0 ; j < tag.size(); j++ )
        cart_vect[i] += tag[j]*lat_vecs[j][i];
    return cart_vect;
};

/**
 * @brief Unit-cell volume.
 * @param uc 3x3 lattice vectors
 * @return volume in cubic Angstrom
 */
double volume( const array < array<double,3> , 3 >&  uc )
{
    return  std::fabs(
            ( uc[0][1]*uc[1][2]-uc[1][1]*uc[0][2])*uc[2][0]+
            ( uc[0][2]*uc[1][0]-uc[0][0]*uc[1][2])*uc[2][1]+
            ( uc[0][0]*uc[1][1]-uc[0][1]*uc[1][0])*uc[2][2]);
}

/**
 * @brief Tight-binding model built from Wannier90 output.
 *
 * Holds the unit cell, orbital positions, and the primitive-cell Hamiltonian as a
 * hopping_list. Provides builder methods for common operators used in KPM/Chebyshev
 * transport calculations.
 */
class tbmodel
{
    public:
    typedef tuple<string, array<double, 3> > orbPos_t;            ///< (label, Cartesian position)
    typedef array < array<double,3> , 3 > unitCell_t;             ///< 3x3 lattice vectors


    /**
     * @brief Read the unit-cell lattice vectors from `.uc`.
     * @param inputfile path to `.uc`
     */
    inline void readUnitCell(const string inputfile)
    {
        lat_vecs= read_unit_cell_file(inputfile);

    }

    /**
     * @brief Read orbital positions from `.xyz`.
     * @param inputfile path to `.xyz`
     */
    inline void readOrbitalPositions(const string inputfile)
    {
        orbPos_list= read_xyz_file(inputfile);
    };

    /**
     * @brief Read the Hamiltonian from `_hr.dat`.
     * @param inputfile path to `_hr.dat`
     */
    inline void readWannierModel(const string inputfile)
    {
        hl = create_hopping_list(read_wannier_file(inputfile));
    };

    /**
     * @brief Apply the Wigner-Seitz minimum-image correction from `_wsvec.dat`.
     * @param inputfile path to `_wsvec.dat`
     *
     * No-op if the file is absent.
     */
    inline void applyWsvec(const string inputfile)
    {
        hl = apply_wsvec(hl, read_wsvec(inputfile));
    };

    /**
     * @brief Load an arbitrary operator in `_hr.dat` format.
     * @param inputfile path to operator file
     * @return hopping_list representation
     *
     * Reuses the same parse + ndegen-normalize chain as the Hamiltonian. Enables
     * hand-built tight-binding models and externally generated operators.
     */
    inline hopping_list readOperatorModel(const string inputfile)
    {
        return create_hopping_list(read_wannier_file(inputfile));
    };

    /**
     * @brief Read the Wannier90 position matrix `_r.dat` (the Berry connection).
     * @param inputfile path to `<seed>_r.dat` (Wannier90 `write_rmn` output)
     * @return array of three hopping_lists, the x/y/z components of
     *         \f$A_{a,ij}(R)=\langle 0i|r_a|Rj\rangle\f$ (Angstrom).
     *
     * `_r.dat` format: a comment line, then num_wann, then nrpts, then one line per
     * (R, m, n): `R1 R2 R3 m n  Re(x) Im(x) Re(y) Im(y) Re(z) Im(z)`. Unlike `_hr.dat`
     * there is NO ndegen block in the file -- Wannier90 reuses the Hamiltonian's
     * ndegen. This reader takes the values verbatim (correct for the common
     * ndegen=1 case and for any file already on H's normalization). Pitfall: for a
     * Wigner-Seitz model with ndegen>1 the position matrix must be pre-normalized to
     * the same convention as the `_hr.dat` consumed here (create_hopping_list divides
     * H by ndegen), or the covariant velocity will be mis-scaled.
     */
    inline array<hopping_list,3> readPositionMatrix(const string inputfile)
    {
        ifstream f(inputfile.c_str());
        if( !f.good() )
            throw std::runtime_error("readPositionMatrix: cannot open '" + inputfile + "'");
        string line;
        std::getline(f, line);                                  // comment line
        int nw=0, nrpts=0;
        { std::getline(f,line); stringstream ss(line); ss>>nw; }
        { std::getline(f,line); stringstream ss(line); ss>>nrpts; }
        if( nw <= 0 )
            throw std::runtime_error("readPositionMatrix: bad num_wann in '" + inputfile + "'");
        array<hopping_list,3> A;
        for( int a=0;a<3;a++ ){ A[a].num_wann = nw; A[a].cellSizes = {1,1,1}; }
        while( std::getline(f, line) )
        {
            stringstream ss(line);
            hopping_list::cellID_t R; for(auto& x:R) ss>>x;
            int m,n; ss>>m>>n; m-=1; n-=1;                      // 1-based -> 0-based
            double re,im;
            for( int a=0;a<3;a++ )
            {
                if( !(ss>>re>>im) ) break;
                hopping_list::value_t v(re,im);
                if( sqrt(norm(v)) > numeric_limits<double>::epsilon() )
                    A[a].hoppings.push_back( hopping_list::hopping_t(R, v, hopping_list::edge_t({m,n})) );
            }
        }
        if( nw != hl.WannierBasisSize() && hl.WannierBasisSize() > 0 )
            throw std::runtime_error("readPositionMatrix: num_wann in _r.dat (" + std::to_string(nw) +
                ") != Hamiltonian num_wann (" + std::to_string(hl.WannierBasisSize()) + ")");
        return A;
    };

    /**
     * @brief Velocity (current) operator along a Cartesian axis.
     * @param dir  direction index: 0=x, 1=y, 2=z
     * @param bare if true, use only the cell displacement R.lat (the gradient
     *             \f$\partial_a H\f$); if false (default), also add the intra-cell
     *             Wannier-centre difference \f$\Delta r_{ij}=r_j-r_i\f$.
     * @return hopping_list representation of \f$v_{dir} = -i\,(\text{displacement})_{dir}\,H\f$
     *
     * This is the velocity ladder's first two rungs (the third, covariant, is
     * createCovariantVelocity_list):
     *   - `bare=true`  : \f$v_a(R) = -i\,(R\cdot\mathrm{lat})_a\,H(R)\f$ — the pure
     *     Bloch-phase gradient, no orbital positions. Cheapest, topology-only.
     *   - `bare=false` : \f$v_a(R) = -i\,(R\cdot\mathrm{lat}+\Delta r_{ij})_a\,H(R)\f$ —
     *     adds the DIAGONAL Berry connection (the Wannier centres from `.xyz`). This
     *     is the historical default and equals \f$v_{\text{bare}}-i[H,A_{\text{diag}}]\f$
     *     because \f$[H,A_{\text{diag}}]_{ij}=(r_j-r_i)H_{ij}\f$.
     *
     * Default (`bare=false`) is byte-identical to the previous one-argument form, so
     * every existing VX/VY golden is preserved. Requires `.uc` and `.xyz`.
     * Units: eV*Angstrom.
     */
    hopping_list createHoppingCurrents_list(const int dir, const bool bare=false)
    {
        std::cout<<"Creating the Current matrix J"<<dir<<(bare?" (bare gradient)":"")<<std::endl;
        assert( dir <3 && dir >=0 );
        assert(volume(lat_vecs) > 0 );
        assert_equal( (int)orbPos_list.size(), hl.WannierBasisSize());

        hopping_list chl = this->hl ;
        for( auto& elem: chl.hoppings )
        {
            auto  tag   = get<0>(elem);
            auto  edge  = get<2>(elem);

            auto displ = tag2cartesian(tag, lat_vecs); //cell displacement R.lat (cartesian)

            if( !bare )                                 //add the Wannier-centre difference r_j - r_i
                for( int i=0; i < displ.size(); i++)
                    displ[i] += get<1>(orbPos_list[edge[1]])[i] - get<1>(orbPos_list[edge[0]])[i];

            //Change the hopping element accordingly
            get<1>(elem) *=  hopping_list::value_t( 0.0, -displ[dir] );
        }
        return chl;
    };

    /**
     * @brief R-space matrix product of two operators: \f$C(R)=\sum_{R'}A(R')\,B(R-R')\f$.
     *
     * Both operands are primitive-cell operators in hopping_list form (the duplicate
     * edges they may carry are summed first). The convolution is the real-space image
     * of the k-space product \f$C(k)=A(k)B(k)\f$, so it is exact at the primitive level
     * (no supercell, no periodic-wrap aliasing). Cost \f$O(nR_A\,nR_B\,n_w^3)\f$ — cheap
     * for primitive cells (nR ~ 1e2, n_w ~ 1e1). Used to form the commutator that turns
     * the bare velocity into the covariant one. The result inherits A's basis size/bounds.
     */
    static hopping_list hopping_product_R(const hopping_list& A, const hopping_list& B)
    {
        typedef hopping_list::value_t cplx;
        typedef hopping_list::cellID_t cell;
        const int nw = A.num_wann;
        auto densify = [nw](const hopping_list& O){
            std::map<cell, std::vector<cplx> > M;
            for( const auto& h : O.hoppings ){
                const auto& R = get<0>(h); const auto& e = get<2>(h);
                auto it = M.find(R);
                if( it == M.end() ) it = M.insert({R, std::vector<cplx>((size_t)nw*nw, cplx(0,0))}).first;
                it->second[(size_t)e[0]*nw + e[1]] += get<1>(h);
            }
            return M;
        };
        const auto MA = densify(A), MB = densify(B);
        std::map<cell, std::vector<cplx> > MC;
        for( const auto& a : MA )
        for( const auto& b : MB )
        {
            const cell R = { a.first[0]+b.first[0], a.first[1]+b.first[1], a.first[2]+b.first[2] };
            auto it = MC.find(R);
            if( it == MC.end() ) it = MC.insert({R, std::vector<cplx>((size_t)nw*nw, cplx(0,0))}).first;
            for( int i=0;i<nw;i++ )
            for( int j=0;j<nw;j++ )
            {
                const cplx aij = a.second[(size_t)i*nw+j];
                if( aij == cplx(0,0) ) continue;
                for( int k=0;k<nw;k++ )
                    it->second[(size_t)i*nw+k] += aij * b.second[(size_t)j*nw+k];
            }
        }
        hopping_list out; out.num_wann = nw; out.cellSizes = A.cellSizes;
        for( const auto& c : MC )
        for( int i=0;i<nw;i++ )
        for( int k=0;k<nw;k++ )
        {
            const cplx v = c.second[(size_t)i*nw+k];
            if( sqrt(norm(v)) > numeric_limits<double>::epsilon() )
                out.hoppings.push_back( hopping_list::hopping_t(c.first, v, hopping_list::edge_t({i,k})) );
        }
        return out;
    }

    /**
     * @brief Covariant velocity \f$v_a = -i(R\cdot\mathrm{lat})_a H - i[H,A_a]\f$.
     * @param dir Cartesian direction (0=x,1=y,2=z)
     * @param A_a position-matrix (Berry connection) operator for this direction,
     *            \f$A_{a,ij}(R)=\langle 0i|r_a|Rj\rangle\f$, from `_r.dat`.
     * @return hopping_list of the covariant velocity (eV*Angstrom)
     *
     * Top rung of the velocity ladder (Wang-Yates-Souza-Vanderbilt, PRB 74, 195118).
     * Built from the BARE gradient plus the FULL Berry-connection commutator, so the
     * inter-Wannier dipoles (the off-diagonal of A_a) enter the interband velocity
     * matrix elements that any Berry-curvature / spin-Hall response is sensitive to.
     * Because [H,A_diag]_{ij}=(r_j-r_i)H_{ij}, passing only the diagonal of A_a here
     * reproduces the berry_connection rung exactly; the off-diagonal part is the new
     * physics. Pitfall: A_a must be on the SAME ndegen-normalized convention as H
     * (the `_r.dat` reader divides by the same ndegen).
     */
    hopping_list createCovariantVelocity_list(const int dir, const hopping_list& A_a)
    {
        std::cout<<"Creating the covariant velocity v"<<dir<<" = -i(R.lat)H - i[H,A]"<<std::endl;
        hopping_list v = this->createHoppingCurrents_list(dir, /*bare=*/true);   // -i(R.lat)_a H
        const hopping_list HA = hopping_product_R(this->hl, A_a);                // H . A_a
        const hopping_list AH = hopping_product_R(A_a, this->hl);                // A_a . H
        // append  -i * [H,A] = -i*(HA - AH)
        for( const auto& h : HA.hoppings )
            v.hoppings.push_back( hopping_list::hopping_t(get<0>(h),
                    hopping_list::value_t(0.0,-1.0)*get<1>(h), get<2>(h)) );
        for( const auto& h : AH.hoppings )
            v.hoppings.push_back( hopping_list::hopping_t(get<0>(h),
                    hopping_list::value_t(0.0,+1.0)*get<1>(h), get<2>(h)) );
        return v;
    }


    /**
     * @brief Density operator: keep onsite terms, zero translations.
     * @return hopping_list representation of the density support
     */
    hopping_list createHoppingDensity_list()
    {
        assert(volume(lat_vecs) > 0 );
        assert(orbPos_list.size()==hl.WannierBasisSize());

        hopping_list chl = this->hl ;
        for( auto& elem: chl.hoppings )
        {
            auto  tag  = get<0>(elem);
            if( tag != hopping_list::cellID_t({0,0,0} ) )//Send to zero  non diagonal elements
                get<1>(elem) *=  0.0 ;
        }

        return chl;
    };

    /**
     * @brief Map orbital index to spin label.
     * @return map index -> +1 for `_s+_`, -1 for `_s-_`, 0 if unmarked
     */
    map<int,int>
    map_id2spin()
    {
        map<int,int> id2spin;
        int id = 0;
        for( auto orb_id : this->orbPos_list )
        {
            auto orb_label = get<0>(orb_id);
            int sz = 0 ;
            if( orb_label.find("_s+_")!=std::string::npos )
                sz= 1;
            if( orb_label.find("_s-_")!=std::string::npos )
                sz=-1;
            if( sz != 0)
                id2spin.insert( {id,sz} );
            id++;
        }
        return id2spin;
    }

    /**
     * @brief Orbital-base label = the orbital label with the sup/sdw spin tag removed.
     */
    static std::string orbital_base_label(std::string s)
    {
        for( const std::string tag : { std::string("_s+_"), std::string("_s-_") } )
        {
            auto p = s.find(tag);
            if( p != std::string::npos ) { s.erase(p, tag.size()); break; }
        }
        return s;
    }

    /**
     * @brief Pair each spinful orbital with its opposite-spin partner (sup/sdw).
     *
     * Spin-operator SUPPORT is built from the spin doubling, not from H's onsite graph.
     * Pairing convention: group by POSITION (primary); orbitals that share a site are
     * disambiguated by the orbital-base label (the label with the _s+_/_s-_ tag removed).
     * Each (position, base-label) group must hold exactly one up (_s+_) and one down
     * (_s-_); otherwise assert -- the model must label its spin channels unambiguously.
     * Spinless orbitals get no partner. The SOC `.spn` path (exact_spin_operator) is
     * untouched and unaffected.
     * @return map orbital index -> opposite-spin partner index (paired orbitals only)
     */
    map<int,int> map_id2partner()
    {
        auto id2spin = this->map_id2spin();
        auto q = [](double x){ return (long long) std::llround(x * 1e6); }; // quantize position
        map< tuple<long long,long long,long long,std::string>, pair<int,int> > grp; // key -> (up_id, down_id)
        int id = 0;
        for( auto& orb : orbPos_list )
        {
            int sz = id2spin.count(id) ? id2spin[id] : 0;
            if( sz != 0 )
            {
                auto& p = get<1>(orb);
                auto key = make_tuple( q(p[0]), q(p[1]), q(p[2]), orbital_base_label(get<0>(orb)) );
                if( grp.find(key) == grp.end() ) grp[key] = std::make_pair(-1,-1);
                int& slot = ( sz > 0 ) ? grp[key].first : grp[key].second;
                // Throw (not assert): assert is a no-op under NDEBUG -> silent-wrong spin
                // operators; a clear exception lets the caller fall back to the .spn route.
                if( slot != -1 )
                    throw std::runtime_error(
                        "map_id2partner: spin pairing ambiguous (>1 orbital with the same "
                        "position+base-label+spin). Use distinct orbital labels, or supply the "
                        "spin operators directly: --op-file SZ <seed>_Sz_hr.dat (from .spn / "
                        "export_operators.py), or --exact-spin.");
                slot = id;
            }
            id++;
        }
        map<int,int> partner;
        for( auto& kv : grp )
        {
            int up = kv.second.first, dn = kv.second.second;
            if( up == -1 || dn == -1 )
                throw std::runtime_error(
                    "map_id2partner: spin pairing incomplete (a position+base-label group lacks "
                    "an up or down partner) -- the label-based spin generator cannot pair this "
                    "model. Provide spin operators via --op-file S{X,Y,Z} <seed>_S{x,y,z}_hr.dat "
                    "(from .spn), or build exact spin with --exact-spin (.spn + _u.mat). A "
                    "tight-binding user can pass any precomputed *_hr.dat the same way.");
            partner[up] = dn;
            partner[dn] = up;
        }
        return partner;
    }

    /**
     * @brief Spin-density operator S_dir (dir = 'x','y','z'), Pauli convention (eig +/-1).
     * @return hopping_list representation
     *
     * The support is built from the SPIN DOUBLING (map_id2partner), NOT from the
     * Hamiltonian's onsite graph. This makes S_alpha a correct Pauli operator even when
     * the onsite block carries same-spin orbital bonds (e.g. graphene's intra-cell
     * sublattice bond, which the old density-support reuse turned into spurious S_z
     * off-diagonals with eigenvalues outside +/-1). Only onsite (R=0) 2x2 blocks per
     * spin pair are emitted:
     *     S_z = diag(+1,-1),  S_x: (up,down)=(down,up)=1,  S_y: (up,down)=-i,(down,up)=+i.
     * Spinless / unpaired orbitals contribute nothing (zero rows/cols).
     */
    hopping_list createHoppingSpinDensity_list(const char dir)
    {
        std::cout<<"Creating the spin Density matrix S"<<dir<<std::endl;
        auto id2spin = this->map_id2spin();
        auto partner = this->map_id2partner();

        hopping_list out = this->hl;   // inherit basis size + cell bounds
        out.hoppings.clear();          // support is rebuilt from the spin doubling
        const hopping_list::cellID_t R0 = {0,0,0};
        typedef hopping_list::value_t cplx;
        typedef hopping_list::edge_t  edge;

        for( auto& kv : partner )
        {
            const int i  = kv.first;       // this orbital
            const int j  = kv.second;      // its opposite-spin partner
            const int si = id2spin[i];     // +1 (up) or -1 (down)
            switch(dir)
            {
                case 'z':  // sigma_z : +/-1 on the diagonal
                    out.hoppings.push_back( make_tuple(R0, cplx((double)si, 0.0), edge({i,i})) );
                    break;
                case 'x':  // sigma_x : 1 on the up<->down off-diagonal (both directions)
                    out.hoppings.push_back( make_tuple(R0, cplx(1.0, 0.0), edge({i,j})) );
                    break;
                case 'y':  // sigma_y : (up,down)=-i, (down,up)=+i  -> value (0,-si)
                    out.hoppings.push_back( make_tuple(R0, cplx(0.0, (double)(-si)), edge({i,j})) );
                    break;
                default:
                    break;
            }
        }
        std::cout<<"Sucess."<<std::endl;
        return out;
    }


    /**
     * @brief Spin-current operator \f$J_{dir} S_{sdir}\f$.
     * @param dir velocity direction index (0=x, 1=y, 2=z)
     * @param sdir spin direction ('x', 'y', 'z')
     * @return hopping_list representation
     *
     * Spin currents are current operators projected through the selected spin Pauli
     * matrix in the same spin-resolved orbital ordering. Returns zero for a spinless
     * model.
     */
    /**
     * @brief Spin-current operator J^{sdir}_{dir} from the velocity and a Pauli spin factor.
     *
     * Multiplies each velocity element v_{ij}(R) by the per-element spin matrix factor.
     * This element-wise form equals the anticommutator J = 1/2{v, sigma} ONLY for the
     * DIAGONAL spin matrix sigma_z: (1/2{v,sigma_z})_{ij} = 1/2 v_{ij}(s_i+s_j), which is
     * v_{ij} s_i for same-spin (s_i=s_j) and ZERO for opposite-spin. The previous code
     * kept the opposite-spin block unchanged instead of zeroing it (a bug for SOC models
     * where v has spin-flip elements); fixed below.
     *
     * Warning: for the OFF-diagonal sigma_x/sigma_y the anticommutator mixes orbital
     * pairs and CANNOT be written as v_{ij} times a per-element factor. The element-wise
     * x/y values below are therefore NOT the true spin current; the authoritative route
     * for x/y is the supercell matrix anticommutator (CSR `--spin-current V S`, which
     * forms 1/2(V*S + S*V)). To-do: build the primitive J for all components as an
     * R-space anticommutator convolution of v(R) and S(R).
     */
    hopping_list createHoppingSpinCurrents_list(const int dir, const char sdir)
    {
        std::cout<<"Creating the spin Current matrix J"<<dir<<"S"<<sdir<<std::endl;
        auto id2spin = this->map_id2spin();
        auto curr = this->createHoppingCurrents_list(dir);
        for( auto& elem: curr.hoppings )
        {
            auto edge  = get<2>(elem);
            auto value =&get<1>(elem);
            auto s1=id2spin[edge[0]], s2= id2spin[edge[1]];
            if( s1!= 0 && s2!= 0 )
            switch(sdir)
            {
                case 'x':                                   // off-diagonal: see Warning above
                *value *= ( s1 + s2 == 0 ) ? hopping_list::value_t(1.0,0.0)
                                           : hopping_list::value_t(0.0,0.0);
                break;

                case 'y':                                   // off-diagonal: see Warning above
                *value *= ( s1 + s2 == 0 ) ? hopping_list::value_t(0.0,(double)s2)
                                           : hopping_list::value_t(0.0,0.0);
                break;

                case 'z':                                   // diagonal: exact anticommutator
                *value *= ( s1 == s2 ) ? hopping_list::value_t((double)s2,0.0)
                                       : hopping_list::value_t(0.0,0.0);  // zero opposite-spin (fix)
                break;

                default:
                    *value *= 0 ;
            }
        }
        std::cout<<"Sucess."<<std::endl;
        return curr;
    }


    int num_orbs;              ///< number of orbitals (redundant with orbPos_list.size())
    unitCell_t lat_vecs;       ///< unit-cell lattice vectors
    vector< orbPos_t > orbPos_list; ///< orbital positions
    hopping_list hl;           ///< primitive-cell Hamiltonian
};

#endif
