#!/bin/bash
# run from project root
cd build
rm -rf _deps/*-subbuild _deps/*-build
rm -f CMakeCache.txt
rm -rf CMakeFiles .cmake hcpwa vendor solver test
echo "Cache cleared. *-src directories preserved."