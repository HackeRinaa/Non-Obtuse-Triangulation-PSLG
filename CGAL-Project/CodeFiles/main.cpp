#include <iostream>
#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include "../HeaderFiles/triangulation_utils.hpp" 

using namespace CGAL;
using namespace std;


// Main function
int main() {
    // Open the input file
    std::ifstream input_file("../data/input.json");
    if (!input_file.is_open()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return 1; // Exit with error
    }

    std::string input_json((std::istreambuf_iterator<char>(input_file)),
                            std::istreambuf_iterator<char>());

    input_file.close();

    if (input_json.empty()) {
        std::cerr << "Error: The input JSON is empty." << std::endl;
        return 1; // Exit with error
    }

    // Parse the JSON input
    nlohmann::json input;
    try {
        input = nlohmann::json::parse(input_json);
    } catch (const nlohmann::json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        return 1; // Exit with error
    }

    std::vector<Point> points;
    std::vector<std::pair<int, int>> constraints;

    parseInput(input, points, constraints);

    CDT cdt;
    for (const auto& p : points) {
        cdt.insert(p);
    }
    if (!constraints.empty()) {
        for (const auto& constraint : constraints) {
            cdt.insert_constraint(points[constraint.first], points[constraint.second]);
        }   
    }

    // Main loop to refine the triangulation by eliminating obtuse angles
    bool improvement = true;
    while (improvement) {
        improvement = false;
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
            if (isObtuse(face)) {
                // First try to resolve the obtuse angle using an edge flip
                CDT::Face_handle face_handle = face;
                bool flipped = applyEdgeFlip(cdt, face_handle); 
                if (!flipped) {
                    // If flipping fails, insert a Steiner point
                    insertSteinerPoint(cdt, face);
                    improvement = true;
                } else {
                    improvement = true;
                }
            }
        }
    }

    
    writeOutput(cdt, points, input["instance_uid"]); 

    return 0;
}
