#include "../HeaderFiles/triangulation_utils.hpp"
#include <fstream>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

using namespace CGAL;
using namespace std;
using json = nlohmann::json;

int main(int argc, char *argv[])
{
    // Ensure correct usage
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input.json> <output.json>" << std::endl;
        return 1;
    }

    // Input and output file paths
    const std::string inputFile = argv[1];
    const std::string outputFile = argv[2];

    // Load input JSON data
    json inputData;
    try
    {
        std::ifstream inFile(inputFile);
        if (!inFile.is_open())
        {
            throw std::runtime_error("Error: Unable to open input file " + inputFile);
        }
        inFile >> inputData;
        std::cout << "Input JSON parsed successfully." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error parsing input JSON: " << e.what() << std::endl;
        return 1;
    }

    // Initialize triangulation and Steiner points
    CDT cdt;
    std::vector<Point> steinerPoints;

    // Perform triangulation based on the method specified in the input
    try
    {
        steinerPoints = performTriangulation(inputData, cdt);
        int obtuseCount = countObtuseTriangles(cdt);
        double maxEdgeLength = evaluateMaxEdgeLength(cdt);

        std::cout << "Triangulation completed successfully." << std::endl;
        std::cout << "Number of Steiner points added: " << steinerPoints.size() << std::endl;
        std::cout << "Obtuse triangles remaining: " << obtuseCount << std::endl;
        std::cout << "Maximum edge length: " << maxEdgeLength << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error during triangulation: " << e.what() << std::endl;
        return 1;
    }

    // Write output results to a JSON file
    try
    {
        writeOutput(inputData, cdt, steinerPoints, outputFile);
        std::cout << "Results written successfully to: " << outputFile << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error writing output JSON: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
