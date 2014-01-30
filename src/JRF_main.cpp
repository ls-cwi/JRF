/*
 * File:   JRF_main.cpp
 * Author: guwek
 *
 * Created on 15 April 2013, 11:57 am
 */

#include <cstdlib>
#include <iostream>
#include <utility>

#include <lemon/list_graph.h>
#include <lemon/arg_parser.h>

#include "Globals.h"
#include "PhylogeneticTree.h"
#include "TreeMetrics.h"

#include <dirent.h>


using namespace std;
using namespace lemon;

/*
 =========================================================================
 TODO
 =========================================================================
 *
 */



int main(int argc, char * const argv[]) {
  // Initialize the argument parser
	ArgParser ap(argc, argv);
  
	string filename_t1, filename_t2, dirname;
  double k = 1.0;
  
	// Add a string option with storage reference for file name
	ap.refOption("t1", "Name of file that contains first tree in Newick format []", filename_t1, false);
  ap.refOption("t2", "Name of file that contains second tree in Newick format []", filename_t2, false);
  ap.refOption("v", "verbosity, 0 = silent, 5 = full [0]", verbosity, false);
	ap.refOption("T", "CPU time limit (s), -1 = no limit [-1]", time_limit, false);
  ap.refOption("d",  "Name of directory with trees in Newick format for all-against-all *.tre comparisons []", dirname, false);
  ap.refOption("k", "Exponent of GRF metric [1.0]", k, false);
	// Perform the parsing process
	// (in case of any error it terminates the program)
	ap.parse();
  
  if (!filename_t1.empty())
  {
    PhylogeneticTree t1, t2;
    t1.parseNewick(filename_t1);
    t2.parseNewick(filename_t2);
    //cout << "t1: " << t1 << endl << endl;
    //cout << "t2: " << t2 << endl;
    
//    cout << "Robinson-Foulds(" << filename_t1 << ", " << filename_t2 << "):\t" << RobinsonFouldsDistance(t1, t2) << endl;
//    cout << "Jaccard-Robinson-Foulds(" << filename_t1 << ", " << filename_t2 << "):\t" << JaccardRobinsonFouldsDistance(t1, t2, k) << endl;
//    std::pair<std::pair<double, double>, double> d = JaccardRobinsonFouldsDistanceCPLEX(t1, t2, k);
//    cout << "Jaccard-Robinson-Foulds_CPLEX(" << filename_t1 << ", " << filename_t2 << "):\t(" << d.first.first << ", " << d.first.second << ") " << d.second << " %" << endl;
//    cout << "|t1|: " << t1.internal_nodes().size() << endl;
//    cout << "|t2|: " << t2.internal_nodes().size() << endl;
    cout <<filename_t1 << "|" << filename_t2 << "\t" << k << "\t" << flush;
    lemon::Timer clk_RF; double d1 = RobinsonFouldsDistance(t1, t2);
    if (verbosity >= 3) cout << d1 << "\t" << (t1.internal_nodes().size() + t2.internal_nodes().size() - d1)/2 << "\t" << clk_RF.userTime() << "\t" << flush;
    lemon::Timer clk_JRF; std::pair<double, int> d2 = JaccardRobinsonFouldsDistance(t1, t2, k);
    if (verbosity >= 3) cout << d2.first << "\t" << d2.second << "\t" << clk_JRF.userTime() << "\t" << flush;
    lemon::Timer clk_JRF_CPLEX; std::pair<std::pair<double, double>, std::pair<double, double> > d3 = JaccardRobinsonFouldsDistanceCPLEX(t1, t2, k);
      if (verbosity >= 3) cout << d3.first.first << "\t" << d3.first.second << "\t" << d3.second.first << "\t" << d3.second.second << "\t" << clk_JRF_CPLEX.userTime() << flush;
    cout << endl;
    
    
  }
  
  if (!dirname.empty())
  {
    if (verbosity >= 5) cout << "[JRF] parsing trees..." << flush;
    DIR*    dir;
    dirent* pdir;
    std::vector<std::string> treefiles;
    
    dir = opendir(dirname.c_str());
    
    while ((pdir = readdir(dir))) {
      string fn = pdir->d_name;
      if (fn.substr(fn.find_last_of(".") + 1) == "tre") treefiles.push_back(fn);
    }
    
    int n = treefiles.size();
    std::vector<PhylogeneticTree*> T(n);
    for (int i = 0; i < n; ++i)
    {
      PhylogeneticTree *t = new PhylogeneticTree();
      t->parseNewick(dirname + "/" + treefiles[i]);
      T[i] = t;
    }
    if (verbosity >= 5) cout << " done (" << treefiles.size() << " trees)." << endl;
        
    if (verbosity >= 5) cout << "[JRF] comparing all-against-all..." << flush;
//    for (int i = 0; i < n; ++i) for (int j = i; j < n; ++j)
//    {
//      //if (verbosity >= 6) cout << "\r" << (i*n+j+1)/(double)(n*(n+1)/2)*100 << " %                      " << flush;
//      cout << treefiles[i] << "|" << treefiles[j] << "\t" << flush;
//      lemon::Timer clk_RF; double d1 = RobinsonFouldsDistance(*T[i], *T[j]);
//      if (verbosity >= 3) cout << d1 << "\t" << clk_RF.userTime() << "\t" << flush;
//      lemon::Timer clk_JRF; double d2 = JaccardRobinsonFouldsDistance(*T[i], *T[j], k);
//      if (verbosity >= 3) cout << d2 << "\t" << clk_JRF.userTime() << "\t" << flush;
//      lemon::Timer clk_JRF_CPLEX; std::pair<std::pair<double, double>, std::pair<double, double> > d3 = JaccardRobinsonFouldsDistanceCPLEX(*T[i], *T[j], k);
//      if (verbosity >= 3) cout << d3.first.first << " " << d3.first.second << " " << d3.second.first << "\t" << clk_JRF_CPLEX.userTime() << endl;
//    }
  }
  
  
  
  return 0;
}
