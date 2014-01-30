#ifndef PHYLO_TREE_H
#define PHYLO_TREE_H

#include <iostream>
#include <lemon/list_graph.h>
#include "Newickform.h"

using namespace lemon;
using namespace std;

class PhylogeneticTree {
  public:
  PhylogeneticTree() : _clade(_t), _label(_t), _dist(_t) {};
  //PhylogeneticTree(const &PhylogeneticTree p) { _t = p.get_graph(); }; // copy constructor
  void parseNewick(std::string filename);
  friend std::ostream& operator <<(std::ostream &os, PhylogeneticTree &t);
  const ListDigraph &get_graph() const { return _t; }
  const ListDigraph::NodeMap<std::string> &getLabels() const { return _label; }
  const ListDigraph::NodeMap<double> &getDistance() const { return _dist; }
//  const list<list<string> > &clades() const { return _clades; }
  const list<ListDigraph::Node> &internal_nodes() const { return _internal_nodes; }
  const list<string> &clade(ListDigraph::Node i) const { return _clade[i]; }
private:
  ListDigraph _t; // the tree structure as directed lemon graph
  ListDigraph::Node _r; // root node
  ListDigraph::NodeMap<std::list<string> > _clade;
  //list<list<string> > _clades;
  //ListDigraph::NodeMap<list<string> > _clade;
  ListDigraph::NodeMap<std::string> _label;
  ListDigraph::NodeMap<double> _dist;
  void _makeTree(newick_node *root, ListDigraph::Node p);
  void _makeClades();
  list<ListDigraph::Node> _internal_nodes;

};

#endif
