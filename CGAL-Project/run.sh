#!/bin/bash
# Remove the 'built' directory if it exists
if [ -d "built" ]; then
  rm -rf built
fi

mkdir built
cd built
cmake ..
make -j4

./triangulation_program ../data/nonConvexParallel/input_non_convex_parallel_aco.json ../data/nonConvexParallel/output_non_convex_parallel_aco.json
./triangulation_program ../data/nonConvexParallel/input_non_convex_parallel_local.json ../data/nonConvexParallel/output_non_convex_parallel_local.json
./triangulation_program ../data/nonConvexParallel/input_non_convex_parallel_sim_an.json ../data/nonConvexParallel/output_non_convex_parallel_sim_an.json
./triangulation_program ../data/nonConvexParallel/input_non_convex_parallel_hybrid.json ../data/nonConvexParallel/output_non_convex_parallel_hybrid.json
