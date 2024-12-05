#include "../HeaderFiles/triangulation_utils.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

using namespace CGAL;
using namespace std;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <input.json> <output.json>" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    json inputData;
    CDT cdt;
    vector<Point> steinerPoints;

    std::cout << "Parsing input JSON..." << std::endl;
    parseInput(inputFile, inputData);

    std::cout << "Performing initial triangulation..." << std::endl;
    json initialResults = performTriangulation(inputData, cdt);

    // Print initial triangulation results
    std::cout << "Initial Triangulation Results:" << std::endl;
    std::cout << initialResults.dump(4) << std::endl; // Pretty print with 4 spaces indentation

    std::cout << "Performing local search for Steiner points..." << std::endl;
    localSearchForSteinerPoints(cdt, steinerPoints, 1000); // Adjust max iterations as needed

    // Extract the updated triangulation results after the local search
    json updatedResults = extractTriangulationResults(cdt);

    // Print updated triangulation results
    std::cout << "Updated Triangulation Results:" << std::endl;
    std::cout << updatedResults.dump(4) << std::endl;

    std::cout << "Writing output to " << outputFile << "..." << std::endl;
    writeOutput(updatedResults, cdt, steinerPoints, outputFile);

    std::cout << "Triangulation and Steiner point resolution completed successfully." << std::endl;
    return 0;
}