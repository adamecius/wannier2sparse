/**
 * @file wannier_parser.hpp
 * @brief Parsers for Wannier90 input files and gauge data.
 *
 * These functions read the raw text files produced by Wannier90 / pw2wannier90:
 * `_hr.dat` (real-space Hamiltonian), `.uc` (lattice vectors), `.xyz` (orbital
 * positions), `.eig` (eigenvalues), `_wsvec.dat` (Wigner-Seitz minimum-image
 * correction), and `.amn` (projection overlaps). Orbital indices are converted to
 * zero-based internal indexing where needed.
 */
#ifndef WANNIER_PARSER
#define WANNIER_PARSER

#include <string>
#include <map>
#include <array>
#include <sstream>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>
#include<limits>
#include <iomanip>
#include <cassert>

using namespace std;

/**
 * @brief Convert a string to int, aborting on failure.
 * @param s input string
 * @return integer value
 */
inline int safe_stoi(const std::string& s ){
    int output;
    try{ output = stoi(s) ;  }
    catch(std::exception const & e)
    {
        cerr<<"error in conversion performed by: " << e.what() <<" in function read_wannier_file"<<endl;
        exit(-1);
    }
    return output;
};

/**
 * @brief Read a Wannier90 `_hr.dat` file.
 * @param wannier_filename path to `_hr.dat`
 * @return tuple of (num_wann, ndegen per R-block, hopping data lines)
 *
 * The hopping lines are returned as raw strings; each line contains
 * `R1 R2 R3 i j Re Im`. ndegen holds the Wigner-Seitz degeneracy weight of each
 * R-point block, in block order.
 */
tuple<int, vector<int>, vector<string> > read_wannier_file(const string wannier_filename);

/**
 * @brief Read a Wannier90-style `.xyz` orbital position file.
 * @param xyz_filename path to `.xyz`
 * @return vector of (orbital_label, Cartesian_position) pairs
 */
vector< tuple<string, array<double, 3> > > read_xyz_file(const string xyz_filename);

/**
 * @brief Read a unit-cell file.
 * @param uc_filename path to `.uc`
 * @return 3x3 array of lattice vectors (rows are lattice vectors)
 */
array< array<double,3> , 3 >  read_unit_cell_file(const string uc_filename);

/**
 * @brief Read a Wannier90 `.eig` file for spectral bounds.
 * @param eig_filename path to `.eig`
 * @param emin output: minimum energy found
 * @param emax output: maximum energy found
 * @return true if the file was read successfully
 *
 * Returns false and leaves emin/emax untouched if the file is missing or empty.
 */
bool read_eig_bounds(const string eig_filename, double& emin, double& emax);

/**
 * @brief One Wigner-Seitz minimum-image correction record.
 *
 * For the hopping (R, i, j), the minimum-image correction replaces it with
 * T.size() copies at R + T, each carrying weight 1/T.size(). Orbital indices are
 * stored zero-based.
 */
struct wsvec_entry
{
    array<int, 3>              R;   ///< original Wigner-Seitz vector
    int                        i, j;///< zero-based orbital indices
    vector< array<int, 3> >    T;   ///< translation shifts for the corrected image(s)
};

/**
 * @brief Parse a Wannier90 `_wsvec.dat` file.
 * @param wsvec_filename path to `_wsvec.dat`
 * @return vector of correction records
 *
 * Format per (R, iw, jw):
 *   "Rx Ry Rz iw jw"  then  "nT"  then nT lines "Tx Ty Tz".
 * Returns empty if the file is missing.
 */
vector<wsvec_entry> read_wsvec(const string wsvec_filename);

#endif
