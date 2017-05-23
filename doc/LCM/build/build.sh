#!/bin/bash -ex

cd "$LCM_DIR"

PACKAGE=$1
BUILD_STRING=$ARCH-$TOOL_CHAIN-$BUILD_TYPE
PREFIX=$PACKAGE-$BUILD_STRING
LOG_FILE="$LCM_DIR/${PREFIX}.log"

ctest -V -S $LCM_DIR/Albany/doc/LCM/build/lcm_build.cmake \
-DSCRIPT_NAME:STRING=`basename $0` \
-DPACKAGE:STRING=$PACKAGE \
-DBUILD_THREADS:STRING=$2 \
| tee $LOG_FILE
