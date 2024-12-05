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
bool isObtuse(const CDT::Face_handle &face)
{
    auto a = face->vertex(0)->point();
    auto b = face->vertex(1)->point();
    auto c = face->vertex(2)->point();

    K::FT ab2 = CGAL::squared_distance(a, b);
    K::FT bc2 = CGAL::squared_distance(b, c);
    K::FT ca2 = CGAL::squared_distance(c, a);

    return ab2 + bc2 < ca2 || ab2 + ca2 < bc2 || bc2 + ca2 < ab2;
}

bool isCloseToExistingPoint(const Point &point, const CDT &cdt, double threshold = 1e-5)
{
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex)
    {
        Point existingPoint = vertex->point();
        if (CGAL::squared_distance(point, existingPoint) < threshold * threshold)
        {
            return true; // Point is close to an existing vertex
        }
    }
    return false; // No close points found
}

json extractTriangulationResults(const CDT &cdt)
{
    json results;
    results["vertices"] = json::array();
    results["edges"] = json::array();
    results["faces"] = json::array();

    // Create a mapping from vertex handles to indices
    std::map<CDT::Vertex_handle, int> vertexIndexMap;
    int index = 0;

    // Extract vertices
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex)
    {
        Point p = vertex->point();
        results["vertices"].push_back({{"x", CGAL::to_double(p.x())}, {"y", CGAL::to_double(p.y())}});
        vertexIndexMap[vertex] = index++;
    }

    // Extract faces and edges
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        json faceData = json::array();
        for (int i = 0; i < 3; ++i)
        {
            Point p = face->vertex(i)->point();
            faceData.push_back({{"x", CGAL::to_double(p.x())}, {"y", CGAL::to_double(p.y())}});
        }
        results["faces"].push_back(faceData);

        // Extract edges as pairs of vertex indices using the mapping
        for (int i = 0; i < 3; ++i)
        {
            int v1_index = vertexIndexMap[face->vertex(i)];
            int v2_index = vertexIndexMap[face->vertex((i + 1) % 3)];
            results["edges"].push_back({v1_index, v2_index});
        }
    }

    return results;
}

int countObtuseTriangles(CDT &cdt, const Point &steinerPoint)
{
    int obtuseCount = 0;

    // Insert the Steiner point temporarily
    CDT::Vertex_handle newVertex = cdt.insert(steinerPoint);

    // Iterate over faces to count obtuse triangles
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        if (isObtuse(face))
        {
            obtuseCount++;
        }
    }

    // Remove the inserted point
    cdt.remove(newVertex);

    return obtuseCount;
}

// Function to perform local search for Steiner points
void localSearchForSteinerPoints(CDT &cdt, vector<Point> &steinerPoints, int maxIterations)
{
    int iterationCount = 0;
    bool improved = true;

    while (improved && iterationCount < maxIterations)
    {
        improved = false;

        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            if (isObtuse(face))
            {
                auto a = face->vertex(0)->point();
                auto b = face->vertex(1)->point();
                auto c = face->vertex(2)->point();

                // Compute circumcenter
                Point circumcenter = CGAL::circumcenter(a, b, c);

                // Check if the circumcenter is a good candidate
                if (isCloseToExistingPoint(circumcenter, cdt))
                {
                    continue; // Skip if it's too close to existing points
                }

                // Count obtuse triangles after adding the circumcenter
                int obtuseCountBefore = countObtuseTriangles(cdt, circumcenter);
                cdt.insert(circumcenter);
                int obtuseCountAfter = countObtuseTriangles(cdt, circumcenter);

                // Compare counts and decide
                if (obtuseCountAfter < obtuseCountBefore)
                {
                    steinerPoints.push_back(circumcenter);
                    improved = true;
                    std::cout << "Inserted Steiner point at: (" << CGAL::to_double(circumcenter.x())
                              << ", " << CGAL::to_double(circumcenter.y()) << ")" << std::endl;
                }
                else
                {
                    cdt.remove(cdt.finite_vertices_begin()); // Remove the circumcenter if it didn't help
                }
            }
        }
        iterationCount++;
    }

    if (iterationCount == maxIterations)
    {
        std::cerr << "Warning: Maximum iteration limit reached in local search for Steiner points." << std::endl;
    }
}

// Function to perform Delaunay triangulation with vertex indexing
json performTriangulation(const json &inputData, CDT &cdt)
{
    json result;

    // Check if points_x and points_y are arrays
    if (!inputData["points_x"].is_array() || !inputData["points_y"].is_array())
    {
        cerr << "Error: points_x and points_y must be arrays." << endl;
        exit(1);
    }

    // Step 1: Insert points into the CDT
    size_t numPoints = inputData["points_x"].size();
    for (size_t i = 0; i < numPoints; ++i)
    {
        double x = inputData["points_x"][i];
        double y = inputData["points_y"][i];
        cdt.insert(Point(x, y));
    }
    std::cout << "Inserted " << cdt.number_of_vertices() << " points into CDT." << std::endl;

    // Step 2: Insert constraints (edges) into the CDT
    if (inputData.contains("additional_constraints"))
    {
        if (!inputData["additional_constraints"].is_array())
        {
            cerr << "Error: additional_constraints must be an array." << endl;
            exit(1);
        }

        const auto &constraints = inputData["additional_constraints"];
        for (const auto &constraint : constraints)
        {
            if (!constraint.is_array() || constraint.size() != 2)
            {
                cerr << "Error: Each constraint must be an array of two indices." << endl;
                exit(1);
            }
            size_t index1 = constraint[0];
            size_t index2 = constraint[1];

            // Ensure indices are within bounds
            if (index1 >= numPoints || index2 >= numPoints)
            {
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

    for (const auto &sp : steinerPoints)
    {
        result["steiner_points_x"].push_back(CGAL::to_double(sp.x()));
        result["steiner_points_y"].push_back(CGAL::to_double(sp.y()));
    }

    result["message"] = "Triangulation completed successfully!";

    return result;
}

Point generateRandomSteinerPoint(CDT &cdt)
{
    // Get the convex hull of the current triangulation (PSLG)
    vector<Point> convexHullPoints;
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex)
    {
        convexHullPoints.push_back(vertex->point());
    }

    // Use the convex hull points to define the boundary for random point generation
    int numPoints = convexHullPoints.size();

    // Ensure there are at least 3 points in the convex hull to form a triangle
    if (numPoints < 3)
    {
        std::cerr << "Error: Not enough points in convex hull to generate random Steiner point." << std::endl;
        exit(1);
    }

    // Randomly select three points from the convex hull to form a triangle
    int idx1 = rand() % numPoints;
    int idx2 = rand() % numPoints;
    int idx3 = rand() % numPoints;

    // Ensure the selected points are distinct
    while (idx1 == idx2 || idx1 == idx3 || idx2 == idx3)
    {
        idx2 = rand() % numPoints;
        idx3 = rand() % numPoints;
    }

    // Get the three vertices of the triangle
    Point p1 = convexHullPoints[idx1];
    Point p2 = convexHullPoints[idx2];
    Point p3 = convexHullPoints[idx3];

    // Generate random barycentric coordinates u and v
    double u = ((double)rand() / RAND_MAX);
    double v = ((double)rand() / RAND_MAX);

    // Ensure the generated point is inside the triangle by using barycentric coordinates
    if (u + v > 1)
    {
        u = 1 - u;
        v = 1 - v;
    }

    // Manually compute the weighted sum (linear combination) of the points
    // using barycentric coordinates
    double w = 1 - u - v;

    // Compute the final Steiner point coordinates
    double x = w * p1.x() + u * p2.x() + v * p3.x();
    double y = w * p1.y() + u * p2.y() + v * p3.y();

    // Return the generated Steiner point
    Point steinerPoint(x, y);
    return steinerPoint;
}

void simulatedAnnealing(CDT &cdt, vector<Point> &steinerPoints, double alpha, double beta, int L)
{
    double temperature = 1.0;
    double coolingRate = 0.99;

    // Compute the initial energy
    int initialObtuseTriangles = countObtuseTriangles(cdt, Point(0, 0));
    int initialSteinerPoints = steinerPoints.size();
    double energy = alpha * initialObtuseTriangles + beta * initialSteinerPoints;

    for (int iteration = 0; iteration < L; ++iteration)
    {
        bool accepted = false;
        for (int trial = 0; trial < 5; ++trial)
        {
            // Select a random Steiner point
            Point candidatePoint = generateRandomSteinerPoint(cdt);

            // Compute the new energy
            int newObtuseTriangles = countObtuseTriangles(cdt, candidatePoint);
            double newEnergy = alpha * newObtuseTriangles + beta * (initialSteinerPoints + 1);

            // Calculate the energy difference
            double deltaEnergy = newEnergy - energy;

            // Metropolis criterion
            if (deltaEnergy < 0 || exp(-deltaEnergy / temperature) > ((double)rand() / RAND_MAX))
            {
                cdt.insert(candidatePoint);              // Insert the new point
                steinerPoints.push_back(candidatePoint); // Add it to the Steiner points list
                energy = newEnergy;                      // Update the energy
                accepted = true;
                break;
            }
        }

        // Update the temperature
        temperature *= coolingRate;

        if (!accepted)
        {
            std::cout << "No better solution found in this iteration." << std::endl;
        }
    }
    std::cout << "Simulated annealing completed." << std::endl;
}

void writeOutput(const json &inputData, const CDT &cdt, const std::vector<Point> &steinerPoints, const std::string &outputFile)
{
    json outputData;

    // Prepare vertices
    outputData["vertices"] = json::array();
    for (auto vertex = cdt.finite_vertices_begin(); vertex != cdt.finite_vertices_end(); ++vertex)
    {
        json vertexData;
        vertexData["x"] = CGAL::to_double(vertex->point().x());
        vertexData["y"] = CGAL::to_double(vertex->point().y());
        outputData["vertices"].push_back(vertexData);
    }

    // Prepare edges
    outputData["edges"] = json::array();
    for (auto edge : cdt.finite_edges())
    {
        json edgeData;
        auto vh1 = edge.first->vertex(edge.second);
        auto vh2 = edge.first->vertex((edge.second + 1) % 3);
        edgeData.push_back(vh1->info().first);
        edgeData.push_back(vh2->info().first);
        outputData["edges"].push_back(edgeData);
    }

    // Prepare faces
    outputData["faces"] = json::array();
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        json faceData = json::array();
        for (int i = 0; i < 3; ++i)
        {
            json vertexData;
            vertexData["x"] = CGAL::to_double(face->vertex(i)->point().x());
            vertexData["y"] = CGAL::to_double(face->vertex(i)->point().y());
            faceData.push_back(vertexData);
        }
        outputData["faces"].push_back(faceData);
    }

    // Prepare Steiner points
    outputData["steiner_points"] = json::array();
    for (const auto &sp : steinerPoints)
    {
        json spData;
        spData["x"] = CGAL::to_double(sp.x());
        spData["y"] = CGAL::to_double(sp.y());
        outputData["steiner_points"].push_back(spData);
    }

    std::cout << "Attempting to write output to " << outputFile << std::endl;
    std::ofstream outFile(outputFile);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open output file " << outputFile << " for writing." << std::endl;
        return;
    }

    outFile << outputData.dump(4);

    if (!outFile)
    {
        std::cerr << "Error: Failed to write data to " << outputFile << "." << std::endl;
    }
    else
    {
        std::cout << "Output successfully written to " << outputFile << "." << std::endl;
    }

    outFile.close();
    if (!outFile.good())
    {
        std::cerr << "Error: Closing the file failed or not all data was flushed." << std::endl;
    }

    std::cout << "File closed successfully. Check " << outputFile << " for output." << std::endl;
}

void parseInput(const string &inputFile, json &inputData)
{
    ifstream inFile(inputFile);
    if (!inFile.is_open())
    {
        cerr << "Error: Unable to open input file " << inputFile << endl;
        exit(1);
    }
    inFile >> inputData;
}