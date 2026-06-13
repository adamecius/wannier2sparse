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

// Returns (num_wann, ndegen[nrpts], hopping_lines). ndegen holds the
// Wigner-Seitz degeneracy weight of each R-point block, in block order.
tuple<int, vector<int>, vector<string> > read_wannier_file(const string wannier_filename);

vector< tuple<string, array<double, 3> > > read_xyz_file(const string xyz_filename);

array< array<double,3> , 3 >  read_unit_cell_file(const string uc_filename);

// Read a Wannier90 seedname.eig file (lines: "band_index kpt_index energy") and
// return the spectral extremes over all entries. Returns false (leaving
// emin/emax untouched) if the file is missing or empty. Optional convenience
// source for the Hamiltonian spectral bounds.
bool read_eig_bounds(const string eig_filename, double& emin, double& emax);

#endif