#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/algorithm.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <nlohmann/json.hpp> // JSON library

using namespace std;
using json = nlohmann::json;

// CGAL type definitions
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef std::pair<int, bool> VertexInfo;
typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_intersections_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> CDT;
typedef CDT::Point Point;

// Function to check if a triangle is obtuse
bool isObtuse(const CDT::Face_handle &face) {
    auto a = face->vertex(0)->point();
    auto b = face->vertex(1)->point();
    auto c = face->vertex(2)->point();

    K::FT ab2 = CGAL::squared_distance(a, b);
    K::FT bc2 = CGAL::squared_distance(b, c);
    K::FT ca2 = CGAL::squared_distance(c, a);

    return ab2 + bc2 < ca2 || ab2 + ca2 < bc2 || bc2 + ca2 < ab2;
}

bool isCloseToExistingPoint(const Point &point, const CDT &cdt, double threshold = 1e-5) {
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        Point existingPoint = vertex->point();
        if (CGAL::squared_distance(point, existingPoint) < threshold * threshold) {
            return true; // Point is close to an existing vertex
        }
    }
    return false; // No close points found
}

json extractTriangulationResults(const CDT &cdt) {
    json results;
    results["vertices"] = json::array();
    results["edges"] = json::array();
    results["faces"] = json::array();

    // Create a mapping from vertex handles to indices
    std::map<CDT::Vertex_handle, int> vertexIndexMap;
    int index = 0;

    // Extract vertices
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        Point p = vertex->point();
        results["vertices"].push_back({{"x", CGAL::to_double(p.x())}, {"y", CGAL::to_double(p.y())}});
        vertexIndexMap[vertex] = index++;
    }

    // Extract faces and edges
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        json faceData = json::array();
        for (int i = 0; i < 3; ++i) {
            Point p = face->vertex(i)->point();
            faceData.push_back({{"x", CGAL::to_double(p.x())}, {"y", CGAL::to_double(p.y())}});
        }
        results["faces"].push_back(faceData);

        // Extract edges as pairs of vertex indices using the mapping
        for (int i = 0; i < 3; ++i) {
            int v1_index = vertexIndexMap[face->vertex(i)];
            int v2_index = vertexIndexMap[face->vertex((i + 1) % 3)];
            results["edges"].push_back({v1_index, v2_index});
        }
    }

    return results;
}



int countObtuseTriangles(CDT &cdt, const Point &steinerPoint) {
    int obtuseCount = 0;

    // Insert the Steiner point temporarily
    CDT::Vertex_handle newVertex = cdt.insert(steinerPoint);
    
    // Iterate over faces to count obtuse triangles
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        if (isObtuse(face)) {
            obtuseCount++;
        }
    }

    // Remove the inserted point
    cdt.remove(newVertex);
    
    return obtuseCount;
}

// Function to perform local search for Steiner points
void localSearchForSteinerPoints(CDT &cdt, vector<Point> &steinerPoints, int maxIterations) {
    int iterationCount = 0;
    bool improved = true;

    while (improved && iterationCount < maxIterations) {
        improved = false;

        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
            if (isObtuse(face)) {
                auto a = face->vertex(0)->point();
                auto b = face->vertex(1)->point();
                auto c = face->vertex(2)->point();

                // Compute circumcenter
                Point circumcenter = CGAL::circumcenter(a, b, c);

                // Check if the circumcenter is a good candidate
                if (isCloseToExistingPoint(circumcenter, cdt)) {
                    continue; // Skip if it's too close to existing points
                }

                // Count obtuse triangles after adding the circumcenter
                int obtuseCountBefore = countObtuseTriangles(cdt, circumcenter);
                cdt.insert(circumcenter);
                int obtuseCountAfter = countObtuseTriangles(cdt, circumcenter);

                // Compare counts and decide
                if (obtuseCountAfter < obtuseCountBefore) {
                    steinerPoints.push_back(circumcenter);
                    improved = true;
                    std::cout << "Inserted Steiner point at: (" << CGAL::to_double(circumcenter.x()) 
                              << ", " << CGAL::to_double(circumcenter.y()) << ")" << std::endl;
                } else {
                    cdt.remove(cdt.finite_vertices_begin()); // Remove the circumcenter if it didn't help
                }
            }
        }
        iterationCount++;
    }

    if (iterationCount == maxIterations) {
        std::cerr << "Warning: Maximum iteration limit reached in local search for Steiner points." << std::endl;
    }
}


// Function to perform Delaunay triangulation with vertex indexing
json performTriangulation(const json &inputData, CDT &cdt) {
    json result;

    // Check if points_x and points_y are arrays
    if (!inputData["points_x"].is_array() || !inputData["points_y"].is_array()) {
        cerr << "Error: points_x and points_y must be arrays." << endl;
        exit(1);
    }

    // Step 1: Insert points into the CDT
    size_t numPoints = inputData["points_x"].size();
    for (size_t i = 0; i < numPoints; ++i) {
        double x = inputData["points_x"][i];
        double y = inputData["points_y"][i];
        cdt.insert(Point(x, y));
    }
    std::cout << "Inserted " << cdt.number_of_vertices() << " points into CDT." << std::endl;

    // Step 2: Insert constraints (edges) into the CDT
    if (inputData.contains("additional_constraints")) {
        if (!inputData["additional_constraints"].is_array()) {
            cerr << "Error: additional_constraints must be an array." << endl;
            exit(1);
        }

        const auto &constraints = inputData["additional_constraints"];
        for (const auto &constraint : constraints) {
            if (!constraint.is_array() || constraint.size() != 2) {
                cerr << "Error: Each constraint must be an array of two indices." << endl;
                exit(1);
            }
            size_t index1 = constraint[0];
            size_t index2 = constraint[1];

            // Ensure indices are within bounds
            if (index1 >= numPoints || index2 >= numPoints) {
                cerr << "Error: Constraint indices out of bounds." << endl;
                exit(1);
            }


            Point p1(inputData["points_x"][index1], inputData["points_y"][index1]);
            Point p2(inputData["points_x"][index2], inputData["points_y"][index2]);
            cdt.insert_constraint(p1, p2);
            std::cout << "Inserted constraint between points: (" << p1.x() << ", " << p1.y() << ") and ("
                      << p2.x() << ", " << p2.y() << ")" << std::endl;
        }
    }

    // Step 3: Perform local search for Steiner points
    vector<Point> steinerPoints;
    localSearchForSteinerPoints(cdt, steinerPoints, 1000); // Adjust max iterations as needed

    // Step 4: Extract triangulation results
    result = extractTriangulationResults(cdt);

    // Step 5: Add metadata to result
    result["num_points"] = cdt.number_of_vertices();
    result["num_edges"] = result["edges"].size(); // Use the edges extracted in extractTriangulationResults
    result["steiner_points_x"] = json::array();
    result["steiner_points_y"] = json::array();

    for (const auto &sp : steinerPoints) {
        result["steiner_points_x"].push_back(CGAL::to_double(sp.x()));
        result["steiner_points_y"].push_back(CGAL::to_double(sp.y()));
    }

    result["message"] = "Triangulation completed successfully!";

    return result;
}

void writeOutput(const json &inputData, const CDT &cdt, const std::vector<Point> &steinerPoints, const std::string &outputFile) {
    json outputData;

    // Prepare vertices
    outputData["vertices"] = json::array();
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex) {
        json vertexData;
        vertexData["x"] = CGAL::to_double(vertex->point().x());
        vertexData["y"] = CGAL::to_double(vertex->point().y());
        outputData["vertices"].push_back(vertexData);
    }

    // Prepare edges
    outputData["edges"] = json::array();
    for (auto edge : cdt.finite_edges()) {
        json edgeData;
        edgeData.push_back(edge.first->vertex(0)->info());
        edgeData.push_back(edge.first->vertex(1)->info());
        outputData["edges"].push_back(edgeData);
    }

    // Prepare faces
    outputData["faces"] = json::array();
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        json faceData = json::array();
        for (int i = 0; i < 3; ++i) {
            json vertexData;
            vertexData["x"] = CGAL::to_double(face->vertex(i)->point().x());
            vertexData["y"] = CGAL::to_double(face->vertex(i)->point().y());
            faceData.push_back(vertexData);
        }
        outputData["faces"].push_back(faceData);
    }

    // Prepare Steiner points
    outputData["steiner_points"] = json::array();
    for (const auto &sp : steinerPoints) {
        json spData;
        spData["x"] = CGAL::to_double(sp.x());
        spData["y"] = CGAL::to_double(sp.y());
        outputData["steiner_points"].push_back(spData);
    }

    // Write to output file
    std::ofstream outFile(outputFile);
    if (!outFile.is_open()) {
        throw std::runtime_error("Error: Could not open output file.");
    }
    outFile << outputData.dump(4); // Pretty print with 4 spaces indentation
    outFile.close();
}

void parseInput(const string &inputFile, json &inputData) { 
    ifstream inFile(inputFile); 
    if (!inFile.is_open()) { 
        cerr << "Error: Unable to open input file " << inputFile << endl; 
        exit(1); 
    } 
    inFile >> inputData; 
}