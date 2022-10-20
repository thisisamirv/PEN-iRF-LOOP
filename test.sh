#!/bin/bash

# Test PEN-iRF-LOOP with example data
echo "Testing PEN-iRF-LOOP:"
PEN_iRF_LOOP_DIR=$PWD
Rscript $PEN_iRF_LOOP_DIR/PEN-iRF-LOOP.R -d $PEN_iRF_LOOP_DIR/Example.RData

echo "Test finished!"