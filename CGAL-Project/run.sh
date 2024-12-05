#!/bin/bash
# Remove the 'built' directory if it exists
if [ -d "built" ]; then
  rm -rf built
fi

mkdir built
cd built
cmake ..
make -j4
./triangulation_program ../data/input.json ../data/output.json
