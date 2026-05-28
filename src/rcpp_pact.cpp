// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "io.h"

//' Run PACT analysis
//'
//' Runs the full PACT (Posterior Analysis of Coalescent Trees) pipeline on a
//' set of phylogenetic trees. Results are written to output files.
//'
//' @param trees_file Path to the NEWICK tree file (e.g., from BEAST or Migrate).
//' @param param_file Path to the PACT parameter file.
//' @param output_prefix Full path prefix for output files. The function appends
//'   \code{.stats}, \code{.tips}, \code{.skylines}, \code{.rules}, and
//'   \code{.pairs} to this prefix.
//'
//' @return Invisibly returns \code{NULL}. Output is written to files.
//'
//' @export
// [[Rcpp::export]]
void pact_run_cpp(std::string trees_file, std::string param_file, std::string output_prefix) {
    IO io(trees_file, param_file, output_prefix);
    io.treeManip();
    io.printTree();
    io.printStatistics();
    io.printTips();
    io.printSkylines();
    io.printPairs();
}
