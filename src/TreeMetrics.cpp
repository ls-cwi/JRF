#include "TreeMetrics.h"
#include "Globals.h"
#include <algorithm>


#include <lemon/matching.h>

// ILOG stuff
#include <ilcplex/ilocplex.h>
#include <ilconcert/iloalg.h>

using namespace std;
using namespace lemon;

ILOLAZYCONSTRAINTCALLBACK6(tree_callback, const ListGraph&, B, const ListGraph::EdgeMap<int>&, ID, const ListGraph::NodeMap<ListGraph::NodeMap<bool>* >&, subset, const ListGraph::NodeMap<ListGraph::NodeMap<bool>* >&, superset, const ListGraph::NodeMap<ListGraph::NodeMap<bool>* >&, disjoint, const IloBoolVarArray&, x)
{
  
  IloEnv env = getEnv();
  //const double epsilon = IloCplex::EpInt;  gunnar: EpInt ist falsch, aber vielleicht was anderes
  //cout << "in callback "  << endl;
  try {
    
    int no_added = 0; // no cuts added yet in this iteration
    
    IloNumArray x_vals(env);
    getValues(x_vals, x);
    
    list<ListGraph::Edge> L;
    for (ListGraph::EdgeIt e(B); e != INVALID; ++e) if (x_vals[ID[e]] > 0) L.push_back(e);
    
    for (list<ListGraph::Edge>::const_iterator e1 = L.begin(); e1 != L.end(); ++e1)
      for (list<ListGraph::Edge>::const_iterator e2 = e1; e2 != L.end(); ++e2)
      {
        //
        //    ListGraph::EdgeIt e2(B);
        //    for (ListGraph::EdgeIt e1(B); e1 != INVALID; ++e1) if (x_vals[ID[e1]] > 1E-5)
        //    {
        //      for (ListGraph::EdgeIt e2 = e1; e2 != INVALID; ++e2) if (x_vals[ID[e2]] > 1E-5)
        //      {
        ListGraph::Node i = B.u(*e1), j = B.v(*e1), k = B.u(*e2), l = B.v(*e2);
        if (i == k) continue;
        if (j == l) continue;
        if (x_vals[ID[*e1]] + x_vals[ID[*e2]] < 1.0) continue;
        if ((*subset[i])[k] && (*subset[j])[l]) continue;
        if ((*superset[i])[k] && (*superset[j])[l]) continue;
        if ((*disjoint[i])[k] && (*disjoint[j])[l]) continue;
        add(x[ID[*e1]] + x[ID[*e2]] <= 1);
        ++no_added;
        
      }
    //cout << "callback added " << no_added << " violated tree constraints" << endl;
    return;
  }
  catch (IloException& e) {
    env.out() << e << endl;
    throw -1;
  }
  catch (Exception& e) {
    env.out() << e.what() << endl;
    throw -1;
  }
  catch (...) {
    env.out() << "something strange in separation callback" << endl;
    throw -1;
  }
}


void print (const list<string> &L)
{
  for (list<string>::const_iterator cit = L.begin(); cit != L.end(); ++cit)
    cout << *cit << " ";
}


inline bool equal(const list<string> &L1, const list<string> &L2)
{
  if (L1.size() != L2.size()) return false;
  if (L1 == L2) return true;
  return false;
}

inline double Jaccard_weight(const list<string> &L1, const list<string> &L2, double k)
{
  if (equal(L1, L2)) return 0.0;
  vector<string>::iterator uit, iit;
  //  cout << L1.size() << "\t" << L2.size() << "\t" << endl;
  //  cout << "L1: "; print(L1); cout << endl;
  //  cout << "L2: "; print(L2); cout << endl;
  vector<string> U(L1.size() + L2.size()), I(min(L1.size(), L2.size()));
  uit = set_union(L1.begin(), L1.end(), L2.begin(), L2.end(), U.begin());
  iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
  U.resize(uit - U.begin());
  I.resize(iit - I.begin());
  //  cout << U.size() << "\t" << I.size() << "\t" << flush;
  //U.unique(); I.unique();
  //cout << "U: "; print(U); cout << endl;
  //cout << "I: "; print(I); cout << endl;
  //cout << "U:  "; print(U); cout << endl;
  //cout << "I:  "; print(I); cout << endl;
  //cout << "2 - 2*(I.size()/(double)U.size())" << 2 - 2*(I.size()/(double)U.size()) << endl;
  //cout << "(" << U.size() << " " << I.size() << ") ";
  //if (I.size()/(double)U.size() < .5) return 10000; danger, danger, high voltage!
  double J = I.size()/(double)U.size();
  //  if (J < .5) return 10000;

  return 2 - 2*pow(J, k);
}

inline double Jaccard_RF_weight(const list<string> &L1, const list<string> &L2)
{
  if (L1.size() != L2.size()) return 10000;
  if (L1 == L2) return 0;
  return 10000;
}

inline double playground_Jaccard_RF_weight(const list<string> &L1, const list<string> &L2)
{
  if (equal(L1, L2)) return 0.0;
  vector<string>::iterator uit, iit;
  vector<string> U(L1.size() + L2.size()), I(min(L1.size(), L2.size()));
  uit = set_union(L1.begin(), L1.end(), L2.begin(), L2.end(), U.begin());
  iit = set_intersection(L1.begin(), L1.end(), L2.begin(), L2.end(), I.begin());
  U.resize(uit - U.begin());
  I.resize(iit - I.begin());
  double J = I.size()/(double)U.size();
  if (J < .5) return 10000;
  return 2 - 2*J;
  
}


inline double Jaccard_RF_weight_unmatched() { return 1.0; }

bool inline subset(const list<string> &A, const list<string> &B)
{
  //cout << "A: "; print(A); cout << endl;
  //cout << "B: "; print(B); cout << endl;
  
  vector<string> D(A.size());
  vector<string>::iterator dit = set_difference(A.begin(), A.end(), B.begin(), B.end(), D.begin());
  D.resize(dit - D.begin());
  //cout << "|D|: " << D.size() << endl;
  
  return D.empty();
}

bool inline disjoint(const list<string> &A, const list<string> &B)
{
  //cout << "A: "; print(A); cout << endl;
  //cout << "B: "; print(B); cout << endl;
  vector<string> I(min(A.size(), B.size()));
  vector<string>::iterator iit = set_intersection(A.begin(), A.end(), B.begin(), B.end(), I.begin());
  I.resize(iit - I.begin());
  //cout << "|I|: " << I.size() << endl;
  return I.empty();
}


bool inline conflict(const list<string> &Y1, const list<string> &Y2, const list<string> &Y1p, const list<string> &Y2p)
{
  if (subset(Y1, Y1p) && subset (Y2, Y2p)) return false;
  if (subset(Y1p, Y1) && subset (Y2p, Y2)) return false;
  if (disjoint(Y1, Y1p) && disjoint(Y2, Y2p)) return false;
  
  return true;
}




double RobinsonFouldsDistance(const PhylogeneticTree &t1, const PhylogeneticTree &t2)
{
  int mismatches = 0;
  //cout << "t1.internal_nodes().size(): " << t1.internal_nodes().size() << endl;
  //cout << "t2.internal_nodes().size(): " << t2.internal_nodes().size() << endl;
  
  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i)
  {
    int mismatch = 1;
    for (list<ListDigraph::Node>::const_iterator j = t2.internal_nodes().begin(); j != t2.internal_nodes().end(); ++j)
      if (equal(t1.clade(*i), t2.clade(*j))) { mismatch = 0; break; }
    mismatches += mismatch;
  }
  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i)
  {
    int mismatch = 1;
    for (list<ListDigraph::Node>::const_iterator j = t1.internal_nodes().begin(); j != t1.internal_nodes().end(); ++j)
      if (equal(t2.clade(*i), t1.clade(*j))) { mismatch = 0; break; }
    mismatches += mismatch;
  }
  
  return mismatches;
}

pair<double, int> JaccardRobinsonFouldsDistance(const PhylogeneticTree &t1, const PhylogeneticTree &t2, double k)
{
  // make bipartite graph with nodes for clades
  ListGraph B; ListGraph::EdgeMap<double> w(B);
  ListDigraph::NodeMap<ListGraph::Node> t1_in_B(t1.get_graph()), t2_in_B(t2.get_graph());
  ListGraph::NodeMap<ListDigraph::Node> B_in_t1(B), B_in_t2(B);
  //ListDigraph::NodeMap<ListGraph::Node> dummy_t1(t1.get_graph()), dummy_t2(t2.get_graph());
  
  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i) { ListGraph::Node v = B.addNode(); t1_in_B[*i] = v; B_in_t1[v] = *i; }
  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i) { ListGraph::Node v = B.addNode(); t2_in_B[*i] = v; B_in_t2[v] = *i; }
  // add dummy "non-matching" nodes
  //  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i) dummy_t1[*i] = B.addNode();
  //  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i) dummy_t2[*i] = B.addNode();
  //  // add edges and weights
  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i)
    for (list<ListDigraph::Node>::const_iterator j = t2.internal_nodes().begin(); j != t2.internal_nodes().end(); ++j)
    {
      w[B.addEdge(t1_in_B[*i], t2_in_B[*j])] = 2 - Jaccard_weight(t1.clade(*i), t2.clade(*j), k);//Jaccard_RF_weight(t1.clade(*i), t2.clade(*j));////playground_Jaccard_RF_weight(t1.clade(*i), t2.clade(*j));//
      //cout << "RF_weight: " << Jaccard_RF_weight(t1.clade(*i), t2.clade(*j)) << ", JRF_weight: " << Jaccard_weight(t1.clade(*i), t2.clade(*j)) << endl;
    }
  //  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i) w[B.addEdge(t1_in_B[*i], dummy_t1[*i])] = -1;
  //  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i) w[B.addEdge(t2_in_B[*i], dummy_t2[*i])] = -1;
  if (verbosity >= 5) cout << "[JRF] computing matching (disrespecting tree structures)... " << flush;
  MaxWeightedMatching<ListGraph, ListGraph::EdgeMap<double> > M(B, w);
  M.run();
  if (verbosity >= 5) cout << "done." << endl;
  //  cout << "weight:    " << M.matchingWeight() << endl;
  //  cout << "t1 int:  " << t1.internal_nodes().size() << endl;
  //  cout << "t2 int:  " << t2.internal_nodes().size() << endl;
  //  cout << "matched: " << M.matchingSize() << endl;
  
  // count conflicts
  // go over pairs of matched edges and count conflicts
  list<ListGraph::Edge> M_e; int conflict_cnt = 0;
  // g: important: only take edges with weight > 0. weight of matching does not change, but all the 0 edges can create problems in terms of conflicts
  for (ListGraph::EdgeIt e(B); e != INVALID; ++e) if (M.matching(e) && w[e] > 1E-4)  M_e.push_back(e);
  for (list<ListGraph::Edge>::const_iterator e1 = M_e.begin(); e1 != M_e.end(); ++e1)
    for (list<ListGraph::Edge>::const_iterator e2 = M_e.begin(); e2 != M_e.end(); ++e2) if (B.id(*e1) < B.id(*e2))
    {
      const list<string> Y1  = t1.clade(B_in_t1[B.u(*e1)]);
      const list<string> Y1p = t1.clade(B_in_t1[B.u(*e2)]);
      const list<string> Y2 =  t2.clade(B_in_t2[B.v(*e1)]);
      const list<string> Y2p = t2.clade(B_in_t2[B.v(*e2)]);
      if (conflict(Y1, Y2, Y1p, Y2p)) ++conflict_cnt;
    }
  
  return make_pair(t1.internal_nodes().size() + t2.internal_nodes().size() - M.matchingWeight(), conflict_cnt);
}

pair<pair<double, double>, pair<double, double> > JaccardRobinsonFouldsDistanceCPLEX(const PhylogeneticTree &t1, const PhylogeneticTree &t2, double k)
{
  IloEnv env; // get ILOG environment
  IloModel M(env); // get model
  IloCplex cplex(env); // get cplex
  
  // shut up cplex
  if (verbosity < 6) {
    cplex.setOut(env.getNullStream());
    cplex.setWarning(env.getNullStream());
    cplex.setError(env.getNullStream());
  }
  // set CPU time limit
  cplex.setParam(IloCplex::ClockType, 1);
  if (time_limit != -1)
    cplex.setParam(IloCplex::TiLim, time_limit);
  
  // zero tolerance (default: 1E-5)
  cplex.setParam(IloCplex::EpGap, 0.0);
  
  
  // make bipartite graph with nodes for clades
  ListGraph B; ListGraph::EdgeMap<double> w(B); ListGraph::EdgeMap<int> ID(B); int edge_cnt = 0;
  ListDigraph::NodeMap<ListGraph::Node> t1_in_B(t1.get_graph()), t2_in_B(t2.get_graph());
  //ListDigraph::NodeMap<ListGraph::Node> dummy_t1(t1.get_graph()), dummy_t2(t2.get_graph());
  
  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i) t1_in_B[*i] = B.addNode();
  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i) t2_in_B[*i] = B.addNode();
  // add dummy "non-matching" nodes
  //  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i) dummy_t1[*i] = B.addNode();
  //  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i) dummy_t2[*i] = B.addNode();
  //  // add edges and weights
  
  ListGraph::NodeMap<ListGraph::NodeMap<bool>* > c_subset(B), c_superset(B), c_disjoint(B);
  for (ListGraph::NodeIt i(B); i != INVALID; ++i)
  {
    c_subset[i] = new ListGraph::NodeMap<bool>(B);
    c_superset[i] = new ListGraph::NodeMap<bool>(B);
    c_disjoint[i] = new ListGraph::NodeMap<bool>(B);
  }
  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i)
    for (list<ListDigraph::Node>::const_iterator j = t1.internal_nodes().begin(); j != t1.internal_nodes().end(); ++j) //if (i != j)
    {
      (*c_subset[t1_in_B[*i]])[t1_in_B[*j]] = subset(t1.clade(*i), t1.clade(*j));
      (*c_superset[t1_in_B[*i]])[t1_in_B[*j]] = subset(t1.clade(*j), t1.clade(*i));
      (*c_disjoint[t1_in_B[*i]])[t1_in_B[*j]] = disjoint(t1.clade(*i), t1.clade(*j));
    }
  
  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i)
    for (list<ListDigraph::Node>::const_iterator j = t2.internal_nodes().begin(); j != t2.internal_nodes().end(); ++j) //if (i != j)
    {
      (*c_subset[t2_in_B[*i]])[t2_in_B[*j]] = subset(t2.clade(*i), t2.clade(*j));
      (*c_superset[t2_in_B[*i]])[t2_in_B[*j]] = subset(t2.clade(*j), t2.clade(*i));
      (*c_disjoint[t2_in_B[*i]])[t2_in_B[*j]] = disjoint(t2.clade(*i), t2.clade(*j));
    }
  
  
  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i)
    for (list<ListDigraph::Node>::const_iterator j = t2.internal_nodes().begin(); j != t2.internal_nodes().end(); ++j)
    {
      ListGraph::Edge e = B.addEdge(t1_in_B[*i], t2_in_B[*j]);
      w[e] = 2 - Jaccard_weight(t1.clade(*i), t2.clade(*j), k);
      ID[e] = edge_cnt++;
      
      //source_in_t1[e] = *i; target_in_t2[e] = *j;
      //cout << "RF_weight: " << Jaccard_RF_weight(t1.clade(*i), t2.clade(*j)) << ", JRF_weight: " << Jaccard_weight(t1.clade(*i), t2.clade(*j)) << endl;
    }
  //  for (list<ListDigraph::Node>::const_iterator i = t1.internal_nodes().begin(); i != t1.internal_nodes().end(); ++i) w[B.addEdge(t1_in_B[*i], dummy_t1[*i])] = -1;
  //  for (list<ListDigraph::Node>::const_iterator i = t2.internal_nodes().begin(); i != t2.internal_nodes().end(); ++i) w[B.addEdge(t2_in_B[*i], dummy_t2[*i])] = -1;
  if (verbosity >= 5) cout << "[JRF] computing matching with CPLEX (disrespecting tree structures)... " << flush;
  
  //cout << t1.internal_nodes().size() * t2.internal_nodes().size() << "\t" << countEdges(B) << endl;
  IloBoolVarArray x(env, t1.internal_nodes().size() * t2.internal_nodes().size());
  
  for (ListGraph::EdgeIt e(B); e != INVALID; ++e) {
    std::stringstream var_name;
    var_name << "x_" << B.id(B.u(e)) << "_" << B.id(B.v(e));
    x[ID[e]] = IloBoolVar(env, var_name.str().c_str());
    M.add(x[ID[e]]);
    //cout << g.id(g.source(a)) << " --> " << g.id(g.target(a)) << "\t" << g.id(a) << " " << g.edge_no(g.id(g.source(a)), g.id(g.target(a))) << endl;
  }
  
  // build objective function
  IloExpr obj_expr(env);
  for (ListGraph::EdgeIt e(B); e != INVALID; ++e) obj_expr += x[ID[e]] * w[e];
  
  M.add(IloObjective(env, obj_expr, IloObjective::Maximize));
  
  // add matching constraints
  for (ListGraph::NodeIt i(B); i != INVALID; ++i)
  {
    IloExpr cons(env);
    for (ListGraph::IncEdgeIt e(B, i); e != INVALID; ++e) cons += x[ID[e]];
    M.add(cons <= 1);
  }
  
  //if (verbosity >= 5) cout << "[JRF]   adding treeish constraints... " << flush;
  
  // add treeish constraints
  cplex.use(tree_callback(env, B, ID, c_subset, c_superset, c_disjoint, x));
  //  ListGraph::EdgeIt e2(B);
  //  for (ListGraph::EdgeIt e1(B); e1 != INVALID; ++e1)
  //  {
  //    for (ListGraph::EdgeIt e2 = e1; e2 != INVALID; ++e2) //if (e1 < e2)
  //    {
  //      ListGraph::Node i = B.u(e1), j = B.v(e1), k = B.u(e2), l = B.v(e2);
  //      if (i == k) continue;
  //      if (j == l) continue;
  //      //cout << "?" << flush;
  //
  //      if ((*c_subset[i])[k] && (*c_subset[j])[l]) continue;
  //      if ((*c_superset[i])[k] && (*c_superset[j])[l]) continue;
  //      if ((*c_disjoint[i])[k] && (*c_disjoint[j])[l]) continue;
  //      M.add(x[ID[e1]] + x[ID[e2]] <= 1);
  //      //cout << "." << flush;
  //    }
  //  }
  //
  //if (verbosity >= 5) cout << "done." << endl;
  
  
  
  if (verbosity >= 5) cout << "[JRF]   solving ILP... " << flush;
  
  cplex.extract(M);
  bool optimal = cplex.solve();
  
  //  if (!optimal) {
  //    cout << endl << endl << endl;
  //    cout << cplex.getStatus() << endl;
  //    cout << endl << endl << endl;
  //    cerr << "JRF: Optimization problems. CPLEX status code " << cplex.getStatus() << endl;
  //    exit(-1);
  //  }
  //
  //if (cplex.getStatus() != IloAlgorithm::Optimal) return -1.0;
  double z = -1.0, z_lower = -1.0, gap = 100000000;
  try {
    z = t1.internal_nodes().size() + t2.internal_nodes().size() - cplex.getObjValue();
    z_lower = t1.internal_nodes().size() + t2.internal_nodes().size() - cplex.getBestObjValue();
    gap = cplex.getMIPRelativeGap()*100;
  } catch (IloException &e)
  { cerr << "cplex error " << cplex.getStatus() << endl;
  }
  
  int no_matched = 0;
  IloNumArray x_vals(env);
  cplex.getValues(x_vals, x);
  for (ListGraph::EdgeIt e(B); e != INVALID; ++e) if (x_vals[ID[e]] > .5 && w[e] > 1E-5) no_matched++;
  
  //clean up
  for (ListGraph::NodeIt i(B); i != INVALID; ++i)
  {
    delete c_subset[i];
    delete c_superset[i];
    delete c_disjoint[i];
  }
  
  return make_pair(make_pair(z_lower, z), make_pair(gap, no_matched));
}


