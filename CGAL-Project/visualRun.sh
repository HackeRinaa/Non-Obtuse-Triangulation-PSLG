#!/bin/bash
# Remove the 'built' directory if it exists
if [ -d "built" ]; then
  rm -rf built
fi

mkdir built
cd built
cmake ..
make

# Run only the visualization executable
./triangulation_visualization
