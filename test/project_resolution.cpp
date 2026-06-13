#include <cassert>
#include <fstream>
#include <string>
#include <vector>
#include <sys/stat.h>
#include <iostream>
#include "w2sp_arguments.hpp"
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
using namespace std;

// parse() takes (argc, argv); build a throwaway argv from strings.
static W2SP_arguments::Status run_parse(W2SP_arguments& a, vector<string> tokens)
{
    vector<char*> argv;
    for (auto& s : tokens) argv.push_back(const_cast<char*>(s.c_str()));
    return a.parse((int)argv.size(), argv.data());
}

int main()
{
    // Default: input prefix is the positional LABEL (positional CLI unchanged).
    {
        W2SP_arguments a;
        assert(run_parse(a, {"prog", "graphene", "2", "2", "1"}) == W2SP_arguments::PROCEED);
        assert(a.label == "graphene");
        assert(a.input_prefix() == "graphene");
    }
    // --seed overrides the input seedname (output still keyed on LABEL).
    {
        W2SP_arguments a;
        assert(run_parse(a, {"prog", "OUT", "1", "1", "1", "--seed", "model"}) == W2SP_arguments::PROCEED);
        assert(a.label == "OUT");
        assert(a.input_prefix() == "model");
    }
    // --project locates inputs in a directory; combines with --seed.
    {
        W2SP_arguments a;
        assert(run_parse(a, {"prog", "OUT", "1", "1", "1", "-p", "dir", "--seed", "m"}) == W2SP_arguments::PROCEED);
        assert(a.input_prefix() == "dir/m");
    }
    {
        W2SP_arguments a;
        assert(run_parse(a, {"prog", "vse", "1", "1", "1", "--project", "samples"}) == W2SP_arguments::PROCEED);
        assert(a.input_prefix() == "samples/vse");
    }

    // Functional: resolve and load a model placed in a subdirectory.
    {
        mkdir("proj", 0755);
        ofstream("proj/m_hr.dat")
            << "chain\n1\n3\n1 1 1\n"
               "-1 0 0 1 1 0.5 0.0\n"
               " 0 0 0 1 1 0.0 0.0\n"
               " 1 0 0 1 1 0.5 0.0\n";

        W2SP_arguments a;
        assert(run_parse(a, {"prog", "m", "4", "1", "1", "--project", "proj"}) == W2SP_arguments::PROCEED);
        assert(a.input_prefix() == "proj/m");

        hopping_list hl = create_hopping_list(read_wannier_file(a.input_prefix() + "_hr.dat"));
        assert(hl.hoppings.size() == 2);   // two nonzero NN hoppings
    }

    std::cout << "PROJECT RESOLUTION TEST PASSED" << std::endl;
    return 0;
}
