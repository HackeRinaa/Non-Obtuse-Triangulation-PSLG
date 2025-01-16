#!/bin/bash
# Remove the 'built' directory if it exists
if [ -d "built" ]; then
  rm -rf built
fi

mkdir built
cd built
cmake ..
make -j4
./triangulation_program ../data/input_aco.json ../data/output_aco.json
./triangulation_program ../data/input_local.json ../data/output_local.json
./triangulation_program ../data/input_hubrid.json ../data/output_hybrid.json