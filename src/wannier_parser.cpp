#include "wannier_parser.hpp"

tuple<int, vector<string> > read_wannier_file(const string wannier_filename)
{
    int num_wann = 0;   // the number of Wannier functions
    int nrpts = 0;      //number of Wigner-Seitz grid-points
    vector< string > wz_gpoints; 
    vector< string > hopping_list; 

    ifstream input_file(wannier_filename.c_str());
    string line; 
    for( int counter = 0; getline(input_file, line); counter++){
        switch( counter ){
            case 0 : continue;  //ignore data or comments
            
            case 1 : num_wann = safe_stoi(line  ); continue;
            
            case 2 : nrpts = safe_stoi(line); continue;
            
            case 3 :{   //Read all Wigner-Seitz grid-points
                const int nrpts_lines = nrpts/15;   //the format impose 15 grid-points per line.  
                wz_gpoints.push_back(line);         //Note the first line is always read so we are always reading nrpts_lines+1
                for(int n=0; n < nrpts_lines ; n++)
                {
                    getline(input_file, line);
                    wz_gpoints.push_back(line);
                    counter++;
                }
               continue;
            }
            default:
                hopping_list.push_back( line );
        }
    }
    input_file.close();        


return tuple<int, vector<string> >(num_wann,hopping_list); 
};


