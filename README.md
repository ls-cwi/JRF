JRF
===

implementation of generalized Robinson-Foulds metric for comparing trees

Compilation instructions
========================

Dependencies
------------

* LEMON 1.3
* ILOG CPLEX (>= 12.0)

Compiling
---------

Get yoshiko from github:

    git clone <HTTPS clone URL (see on the right side of this page)>


First, LEMON 1.3 needs to be installed:

    wget http://lemon.cs.elte.hu/pub/sources/lemon-1.3.tar.gz
    tar xvzf lemon-1.3.tar.gz
    cd lemon-1.3
    cmake -DCMAKE_INSTALL_PREFIX=~/lemon
    make install
    
Note: On Mac OS 10.9, comment out the following two lines and add the code below at line 162 in `CMakeLists.txt` before `make install`

    #ADD_SUBDIRECTORY(demo) 
    #ADD_SUBDIRECTORY(tools)
    
    if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
      set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
    endif()

You can remove the LEMON sources now, i.e., `rm -rf lemon-1.3`. 

CPLEX is a commercial product owned by IBM. For academic purposed it can be obtained at no charge via IBM's Academic Initiative programme:

  http://www-03.ibm.com/ibm/university/academic/pub/page/membership

Next, JRF can be compiled:

    mkdir build
    cd build
    cmake ..
    make

In case auto-detection of LEMON or CPLEX fails, do

    cmake \
    -DLIBLEMON_ROOT=~/lemon \
    -DCPLEX_INC_DIR=~/ILOG/cplex/include/ \
    -DCPLEX_LIB_DIR=~/ILOG/cplex/lib/x86-64_osx/static_pic \
    -DCONCERT_LIB_DIR=~/ILOG/concert/lib/x86-64_osx/static_pic \
    -DCONCERT_INC_DIR=~/ILOG/concert/include/ ..

Running JRF
=============

To run JRF on two small toy trees:

    ./JRF -t1 ../data/t1.txt -t2 ../data/t2.txt -v 4

This will create output similar to

   ../data/t1.txt|../data/t2.txt	1	16	0	0	3.99127	7	0	4.65794	4.65794	0	7	0

The numbers have the following meaning (in order of occurence):
   1:       value of k 
   16:      Robinson-Foulds distance (RF)
   0:       #internal nodes in t1 + #internal nodes in t2 -RF
   0:       time to compute RF
   3.99127: cost of minimum cost bipartite matching (MCBM)
   7:       number of conflicts induced by MCBM
   0:       time to compute MCBM
   4.65794: lower bound on minimum cost arboreal bipartite matching (MCABM)
   4.65794: upper bound on MCABM
   0:       optimality gap
   7:       number of clade pairs matched by MCABM 
   0:       time to compute MCABM 

Get a list of options:

    ./JRF -h
