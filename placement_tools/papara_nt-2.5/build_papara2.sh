#!/bin/sh

# minimal build script for papara 2.0.
# please use the supplied cmake files for anything more than just building papara 2.0.



g++ -o papara -O3 -msse4a -std=c++11 -I. -I ivy_mike/src/ -I ublasJama-1.0.2.3 papara.cpp pvec.cpp pars_align_seq.cpp pars_align_gapp_seq.cpp parsimony.cpp sequence_model.cpp papara2_main.cpp blast_partassign.cpp align_utils.cpp ivy_mike/src/time.cpp ivy_mike/src/tree_parser.cpp ivy_mike/src/getopt.cpp ivy_mike/src/demangle.cpp ivy_mike/src/multiple_alignment.cpp ublasJama-1.0.2.3/EigenvalueDecomposition.cpp -lpthread

#-I/usr/include/boost141/

# if you are using centos 5 with boost 1.41 from the EPEL repository add the following to the command line:
# -I/usr/include/boost141/

# Otherwise if you get compile errors related to boost, download+unpack the latest boost version from boost.org and add -I<boost dir> to the commandline.

