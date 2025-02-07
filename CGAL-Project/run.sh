#!/bin/bash
# Remove the 'built' directory if it exists
if [ -d "built" ]; then
  rm -rf built
fi

mkdir built
cd built
cmake ..
make -j4

./triangulation_program ../data/convexClosedConstraints/input_convex_closed_constraints_hybrid.json ../data/convexClosedConstraints/output_convex_closed_constraints_hybrid.json
