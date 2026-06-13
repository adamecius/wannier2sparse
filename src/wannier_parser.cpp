#include "wannier_parser.hpp"

tuple<int, vector<int>, vector<string> > read_wannier_file(const string wannier_filename)
{
    ifstream input_file(wannier_filename.c_str());
    assert(input_file.is_open());
    input_file.precision( numeric_limits<double>::digits10+2);

    string line;
    getline(input_file, line);                                 // line 0: comment / date
    getline(input_file, line); const int num_wann = safe_stoi(line); assert(num_wann>0);
    getline(input_file, line); const int nrpts    = safe_stoi(line); assert(nrpts>0);

    // Wigner-Seitz degeneracies: nrpts integers, 15 per line. Read token by
    // token until exactly nrpts are collected. This is robust to the case where
    // nrpts is a multiple of 15: a line-based count (nrpts/15 + 1) would consume
    // an extra line and eat the first hopping, desynchronizing everything.
    vector<int> ndegen; ndegen.reserve(nrpts);
    int deg;
    while( (int)ndegen.size() < nrpts && (input_file >> deg) )
        ndegen.push_back(deg);
    assert( (int)ndegen.size() == nrpts );
    getline(input_file, line);                                 // consume rest of the last degeneracy line

    vector< string > hopping_lines;
    while( getline(input_file, line) )
        hopping_lines.push_back(line);

    input_file.close();

return make_tuple(num_wann, ndegen, hopping_lines);
};

vector< tuple<string, array<double, 3> > > read_xyz_file(const string xyz_filename)
{
    typedef tuple<string, array<double, 3> > xyz_elem;
    vector< xyz_elem > xyz_data;
    std::cout<<"Opening File: "<<xyz_filename.c_str()<<std::endl;
    ifstream input_file(xyz_filename.c_str()); assert(input_file.is_open());
    input_file.precision( numeric_limits<double>::digits10+2);
    int num_sites;
    input_file>>num_sites; assert(num_sites>0);

    std::string label;
    array<double, 3>  pos;
    for( int i = 0; i < num_sites; i++)
    {
       input_file>>label>>pos[0]>>pos[1]>>pos[2];
       xyz_data.push_back(xyz_elem(label,pos) );
    } 
    return xyz_data;
};

array< array<double,3> , 3 >  read_unit_cell_file(const string uc_filename)
{
    constexpr int DIM = 3;
    ifstream input_file(uc_filename.c_str()); assert(input_file.is_open());
    input_file.precision( numeric_limits<double>::digits10+2);
    array< array<double,DIM> , DIM > unit_cell;
    for( auto & lat_vec : unit_cell )
    for( auto & li : lat_vec )
       input_file>>li;
    return unit_cell;

};


