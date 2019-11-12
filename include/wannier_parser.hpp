#ifndef WANNIER_PARSER
#define WANNIER_PARSER

#include <string>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

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

tuple<int, vector<string> > read_wannier_file(const string wannier_filename);



#endif