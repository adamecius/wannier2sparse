// Phase B3: parse a Wannier90 .win file into WannierProvenance for the manifest.
// Verifies the documented keys (num_wann/num_bands/mp_grid/exclude_bands/the
// projections block/disentanglement windows/use_ws_distance) parse correctly,
// that comments and `=`/`:` separators and Fortran logicals are handled, and that
// a file missing optional keys degrades gracefully (no throw, flags stay false).
#include <cassert>
#include <fstream>
#include <string>
#include <iostream>
#include "win_parser.hpp"

using namespace std;

int main()
{
    // A representative .win with comments, both separators, and a projections block.
    {
        ofstream("full.win") <<
            "! graphene wannierisation\n"
            "num_wann  = 4\n"
            "num_bands : 8\n"
            "exclude_bands = 1-2,5\n"
            "mp_grid = 6 6 1   ! Monkhorst-Pack\n"
            "dis_win_min = -5.0\n"
            "dis_win_max = 12.0\n"
            "dis_froz_max = 2.0\n"
            "use_ws_distance = .true.\n"
            "\n"
            "begin projections\n"
            "C : pz\n"
            "f=0.0,0.0,0.0 : s ; p\n"
            "end projections\n"
            "\n"
            "begin kpoint_path\n"
            "G 0.0 0.0 0.0   M 0.5 0.0 0.0\n"
            "M 0.5 0.0 0.0   K 0.3333333 0.3333333 0.0\n"
            "K 0.3333333 0.3333333 0.0   G 0.0 0.0 0.0\n"
            "end kpoint_path\n"
            "\n"
            "begin unit_cell_cart\n"
            "2.46 0.0 0.0\n"
            "end unit_cell_cart\n";

        WannierProvenance w = parse_win_provenance("full.win");
        assert(w.present);
        assert(w.source_file == "full.win");
        assert(w.has_num_wann && w.num_wann == 4);
        assert(w.has_num_bands && w.num_bands == 8);
        assert(w.exclude_bands.size() == 3 &&
               w.exclude_bands[0] == 1 && w.exclude_bands[1] == 2 && w.exclude_bands[2] == 5);
        assert(w.has_mp_grid && w.mp_grid[0] == 6 && w.mp_grid[1] == 6 && w.mp_grid[2] == 1);
        assert(w.has_dis_win && w.dis_win_min == -5.0 && w.dis_win_max == 12.0);
        assert(w.has_dis_froz && w.dis_froz_max == 2.0);
        assert(w.use_ws_distance == true);
        assert(w.dis_enabled);   // num_bands > num_wann and windows present
        // Projections block: two entries, the unit_cell_cart block is NOT captured.
        assert(w.projections.size() == 2);
        assert(w.projections[0].site == "C" && w.projections[0].orbital == "pz");
        assert(w.projections[1].site == "f=0.0,0.0,0.0");
        assert(w.projections[1].orbital == "s ; p");
        // kpoint_path: three segments G-M, M-K, K-G collapse to the node list G-M-K-G.
        assert(w.kpoint_path.size() == 4);
        assert(w.kpoint_path[0].label == "G" && w.kpoint_path[1].label == "M" &&
               w.kpoint_path[2].label == "K" && w.kpoint_path[3].label == "G");
        assert(w.kpoint_path[1].k[0] == 0.5 && w.kpoint_path[1].k[1] == 0.0);
        assert(w.kpoint_path_source == "win");
    }

    // Minimal .win: only num_wann. Everything else stays unset; no throw.
    {
        ofstream("min.win") << "num_wann = 2\n";
        WannierProvenance w = parse_win_provenance("min.win");
        assert(w.present);
        assert(w.has_num_wann && w.num_wann == 2);
        assert(!w.has_num_bands && !w.has_mp_grid && !w.has_dis_win && !w.has_dis_froz);
        assert(w.exclude_bands.empty() && w.projections.empty());
        assert(w.kpoint_path.empty() && w.kpoint_path_source.empty());
        assert(w.use_ws_distance == false);
        assert(!w.dis_enabled);
    }

    // use_ws_distance false form, and a missing file must throw.
    {
        ofstream("false.win") << "num_wann = 1\nuse_ws_distance = .false.\n";
        WannierProvenance w = parse_win_provenance("false.win");
        assert(w.use_ws_distance == false);

        bool threw = false;
        try { parse_win_provenance("does_not_exist.win"); } catch (const std::exception&) { threw = true; }
        assert(threw && "a missing .win must throw");
    }

    // QE bands.in (K_POINTS crystal_b) fallback k-path parser.
    {
        ofstream("bands.in") <<
            "&control\n  calculation = 'bands'\n/\n"
            "K_POINTS crystal_b\n"
            "3\n"
            "  0.0 0.0 0.0 20 ! G\n"
            "  0.5 0.0 0.0 20 ! M\n"
            "  0.3333333 0.3333333 0.0 1 ! K\n";
        std::vector<KpathNode> p = parse_qe_bands_kpath("bands.in");
        assert(p.size() == 3);
        assert(p[0].label == "G" && p[1].label == "M" && p[2].label == "K");
        assert(p[1].k[0] == 0.5);

        bool threw = false;
        try { parse_qe_bands_kpath("does_not_exist.in"); } catch (const std::exception&) { threw = true; }
        assert(threw && "a missing QE bands file must throw");
    }

    cout << "WIN PROVENANCE TEST PASSED" << endl;
    return 0;
}
