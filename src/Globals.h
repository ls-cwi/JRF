//
//  Globals.h
//  JRF
//
//  Created by Gunnar W. Klau on 15-4-13.
//
//

#ifndef __JRF__Globals__
#define __JRF__Globals__

#include <limits>
#include <cstdlib>

#include <lemon/time_measure.h>


extern int time_limit;
extern lemon::Timer clk;
extern int verbosity;
extern int no_threads;
extern double eps;

#endif /* defined(__JRF__Globals__) */
