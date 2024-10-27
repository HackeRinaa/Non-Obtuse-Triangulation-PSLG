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

// Main function
int main()
{
    // Open the input file
    std::ifstream input_file("../data/input.json");
    if (!input_file.is_open())
    {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1; // Exit with error
    }

    std::string input_json((std::istreambuf_iterator<char>(input_file)),
                           std::istreambuf_iterator<char>());
    input_file.close();

    if (input_json.empty())
    {
        std::cerr << "Error: The input JSON is empty." << std::endl;
        return 1; // Exit with error
    }

    // Parse the JSON input
    nlohmann::json input;
    try
    {
        input = nlohmann::json::parse(input_json);
    }
    catch (const nlohmann::json::parse_error &e)
    {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        std::cerr << "Input JSON: " << input_json << std::endl; // Show what was attempted to parse
        return 1;                                               // Exit with error
    }

    std::vector<Point> points;
    std::vector<std::pair<int, int>> constraints;

    try
    {
        parseInput(input, points, constraints);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error while parsing input: " << e.what() << std::endl;
        return 1; // Exit with error
    }

    // Check if points and constraints are populated
    if (points.empty())
    {
        std::cerr << "Error: No points were parsed from the input." << std::endl;
        return 1; // Exit with error
    }

    if (constraints.empty())
    {
        std::cout << "No additional constraints provided." << std::endl;
    }
    else
    {
        std::cout << "Parsed constraints:" << std::endl;
        for (const auto &constraint : constraints)
        {
            std::cout << "Constraint: (" << constraint.first << ", " << constraint.second << ")" << std::endl;
        }
    }

    CDT cdt;
    for (const auto &p : points)
    {
        cdt.insert(p);
    }
    if (!constraints.empty())
    {
        for (const auto &constraint : constraints)
        {
            cdt.insert_constraint(points[constraint.first], points[constraint.second]);
        }
    }
    bool improvement = true;
    while (improvement)
    {
        improvement = false;
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            if (isObtuse(face))
            {
                // First try to resolve the obtuse angle using an edge flip
                CDT::Face_handle face_handle = face;
                bool flipped = applyEdgeFlip(cdt, face_handle);
                if (!flipped)
                {
                    std::cerr << "Inserting Steiner point for face at: " << face->vertex(0)->point() << ", "
                              << face->vertex(1)->point() << ", " << face->vertex(2)->point() << std::endl;
                    insertSteinerPoint(cdt, face);
                    improvement = true;
                }
                else
                {
                    std::cerr << "Flipped an edge for face: " << face->vertex(0)->point() << ", "
                              << face->vertex(1)->point() << ", " << face->vertex(2)->point() << std::endl;
                    improvement = true;
                }
            }
        }

        // Prepare output
        nlohmann::json output;
        writeOutput(cdt, points, output);

        // Output results to console or file
        std::ofstream output_file("../data/output.json");
        if (output_file.is_open())
        {
            output_file << output.dump(4); // Pretty print with 4 spaces
            output_file.close();
            std::cout << "Output written to ../data/output.json" << std::endl;
        }
        else
        {
            std::cerr << "Error: Could not open output file." << std::endl;
            return 1; // Exit with error
        }

        return 0;
    }
}