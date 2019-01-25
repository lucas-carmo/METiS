#!/bin/bash
# Run METiS - REMOVE THIS PART AFTER TESTS ARE DONE
#
#

# Need that to find libhdf5.so.101
LD_LIBRARY_PATH=/home/lucas/anaconda3/lib
export LD_LIBRARY_PATH

./bin/METiS "$1"