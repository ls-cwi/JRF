#ifndef TREE_METRICS_H
#define TREE_METRICS_H

#include <iostream>
#include "PhylogeneticTree.h"

using namespace lemon;

double RobinsonFouldsDistance(const PhylogeneticTree &t1, const PhylogeneticTree &t2);
pair<double, int> JaccardRobinsonFouldsDistance(const PhylogeneticTree &t1, const PhylogeneticTree &t2, double k);
pair<pair<double, double>, pair<double, double> > JaccardRobinsonFouldsDistanceCPLEX(const PhylogeneticTree &t1, const PhylogeneticTree &t2, double k);

#endif
