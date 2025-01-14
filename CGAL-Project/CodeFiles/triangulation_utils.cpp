#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/algorithm.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <cmath>
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

    // Check if any angle is obtuse using the squared lengths
    return ab2 + bc2 < ca2 || ab2 + ca2 < bc2 || bc2 + ca2 < ab2;
}

// Function to check if the triangle has an obtuse angle
bool hasObtuseAngle(const Point &p0, const Point &p1, const Point &p2)
{
    double a2 = CGAL::squared_distance(p1, p2); // distance squared
    double b2 = CGAL::squared_distance(p0, p2);
    double c2 = CGAL::squared_distance(p0, p1);

    return (b2 + c2 < a2 || a2 + c2 < b2 || a2 + b2 < c2); // Check for obtuseness
}

int countObtuseTriangles(const CDT &cdt)
{
    int obtuseCount = 0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        if (cdt.is_infinite(face))
            continue; // Skip infinite faces

        Point p0 = face->vertex(0)->point();
        Point p1 = face->vertex(1)->point();
        Point p2 = face->vertex(2)->point();

        if (hasObtuseAngle(p0, p1, p2))
        {
            obtuseCount++;
        }
    }
    return obtuseCount;
}

Point getMidpoint(const Point &p1, const Point &p2)
{
    return Point((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2);
}

Point getCentroid(const Point &p0, const Point &p1, const Point &p2)
{
    return Point((p0.x() + p1.x() + p2.x()) / 3, (p0.y() + p1.y() + p2.y()) / 3);
}

Point getCircumcenter(const Point &p0, const Point &p1, const Point &p2)
{
    // Use formula to calculate the circumcenter
    double D = 2 * (p0.x() * (p1.y() - p2.y()) + p1.x() * (p2.y() - p0.y()) + p2.x() * (p0.y() - p1.y()));
    double Ux = ((p0.x() * p0.x() + p0.y() * p0.y()) * (p1.y() - p2.y()) + (p1.x() * p1.x() + p1.y() * p1.y()) * (p2.y() - p0.y()) + (p2.x() * p2.x() + p2.y() * p2.y()) * (p0.y() - p1.y())) / D;
    double Uy = ((p0.x() * p0.x() + p0.y() * p0.y()) * (p2.x() - p1.x()) + (p1.x() * p1.x() + p1.y() * p1.y()) * (p0.x() - p2.x()) + (p2.x() * p2.x() + p2.y() * p2.y()) * (p1.x() - p0.x())) / D;

    return Point(Ux, Uy);
}

Point projectVertexOntoEdge(const Point &p0, const Point &p1, const Point &p2)
{
    // Compute the projection of p0 onto the edge p1p2
    double x0 = p0.x(), y0 = p0.y();
    double x1 = p1.x(), y1 = p1.y();
    double x2 = p2.x(), y2 = p2.y();

    double dx = x2 - x1;
    double dy = y2 - y1;

    double t = ((x0 - x1) * dx + (y0 - y1) * dy) / (dx * dx + dy * dy);

    // The projection point
    double px = x1 + t * dx;
    double py = y1 + t * dy;

    return Point(px, py);
}

double sign(const Point &p1, const Point &p2, const Point &p3)
{
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
}

bool isPointInsideTriangle(const Point &p, const Point &p0, const Point &p1, const Point &p2)
{
    double d1 = sign(p, p0, p1);
    double d2 = sign(p, p1, p2);
    double d3 = sign(p, p2, p0);

    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos); // If all signs are the same, point is inside
}

std::vector<Point> generateCandidatePoints(const Point &p1, const Point &p2, const Point &p3)
{
    std::vector<Point> candidates;

    // Circumcenter
    candidates.push_back(CGAL::circumcenter(p1, p2, p3));

    // Midpoints of edges
    candidates.push_back(CGAL::midpoint(p1, p2));
    candidates.push_back(CGAL::midpoint(p2, p3));
    candidates.push_back(CGAL::midpoint(p3, p1));

    // Centroid
    candidates.push_back(CGAL::centroid(p1, p2, p3));

    candidates.push_back(generateRandomPointInTriangle(p1, p2, p3));

    return candidates;
}

// Function to evaluate the triangulation and penalize obtuse triangles
double evaluateTriangulation(const CDT &cdt)
{
    double score = 0.0;
    int totalTriangles = 0;
    int obtuseCount = 0;

    // Count the total number of triangles and the obtuse triangles
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        if (cdt.is_infinite(face))
            continue;

        Point p0 = face->vertex(0)->point();
        Point p1 = face->vertex(1)->point();
        Point p2 = face->vertex(2)->point();

        totalTriangles++;

        if (hasObtuseAngle(p0, p1, p2))
        {
            obtuseCount++;
        }
    }

    // Calculate the score as a fraction of obtuse triangles
    if (totalTriangles > 0)
    {
        score = -static_cast<double>(obtuseCount) / totalTriangles;
    }

    // Debugging output
    cout << "Total triangles: " << totalTriangles << ", Obtuse triangles: " << obtuseCount << ", Score: " << score << endl;

    return score;
}

// Function to perform local search for Steiner points
void localSearchForSteinerPoints(CDT &cdt, std::vector<Point> &steinerPoints, int maxIterations)
{
    for (int iter = 0; iter < maxIterations; ++iter)
    {
        bool improvementMade = false;
        cout << "Iteration: " << iter << endl;

        // Iterate through all finite faces (triangles) in the triangulation
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
        {
            if (cdt.is_infinite(face))
                continue; // Skip infinite faces

            // Get the vertices of the triangle
            Point p0 = face->vertex(0)->point();
            Point p1 = face->vertex(1)->point();
            Point p2 = face->vertex(2)->point();

            // Skip triangles that do not have obtuse angles
            if (!hasObtuseAngle(p0, p1, p2))
                continue;

            // Generate candidate Steiner points
            std::vector<Point> candidatePoints = generateCandidatePoints(p0, p1, p2);
            if (candidatePoints.empty())
            {
                cout << "No candidate points generated for triangle (" << p0 << ", " << p1 << ", " << p2 << ")" << endl;
                continue;
            }

            Point bestPoint;
            int minObtuseCount = std::numeric_limits<int>::max(); // Start with a large number of obtuse triangles

            // Evaluate each candidate point
            for (const auto &candidate : candidatePoints)
            {

                // Insert the candidate point and check if the insertion is successful
                auto vh = cdt.insert(candidate);
                if (vh == nullptr)
                {
                    continue; // Skip this candidate if insertion fails
                }

                // Count the number of obtuse triangles after insertion
                int obtuseCount = countObtuseTriangles(cdt);

                // Check if this candidate results in fewer obtuse triangles
                if (obtuseCount < minObtuseCount)
                {
                    bestPoint = candidate;
                    minObtuseCount = obtuseCount;
                }

                // Remove the candidate point to evaluate the next one
                cdt.remove(vh);
            }

            // If a valid candidate is found that reduces obtuse triangles, insert it
            int currentObtuseCount = countObtuseTriangles(cdt); // Get current obtuse count
            if (minObtuseCount < currentObtuseCount)
            {                                       // Compare with current obtuse count
                cdt.insert(bestPoint);              // Insert the best point found
                steinerPoints.push_back(bestPoint); // Add best point to the list
                cout << "Inserted Steiner point: " << bestPoint << " with obtuse triangles: " << minObtuseCount << endl;
                improvementMade = true;
            }
            else
            {
                cout << "No improvement found for triangle (" << p0 << ", " << p1 << ", " << p2 << ")" << endl;
            }
        }

        // Exit the loop if no improvement was made
        if (!improvementMade)
            break;
    }
}

// Function to extract triangulation results into JSON format
json extractTriangulationResults(const CDT &cdt)
{
    json result;
    result["triangles"] = json::array();

    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        if (cdt.is_infinite(face))
            continue;

        json triangle;
        for (int i = 0; i < 3; ++i)
        {
            Point p = face->vertex(i)->point();
            triangle.push_back({CGAL::to_double(p.x()), CGAL::to_double(p.y())});
        }
        result["triangles"].push_back(triangle);
    }

    return result;
}

// Function to initialize pheromones
std::map<Point, double> initializePheromones(const CDT &cdt)
{
    std::map<Point, double> pheromones;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        if (cdt.is_infinite(face))
            continue;
        Point p0 = face->vertex(0)->point();
        Point p1 = face->vertex(1)->point();
        Point p2 = face->vertex(2)->point();
        pheromones[getCentroid(p0, p1, p2)] = 1.0; // Initialize to 1.0
    }
    return pheromones;
}

// Function to calculate heuristic value
double calculateHeuristic(const CDT::Face_handle &face)
{
    Point p0 = face->vertex(0)->point();
    Point p1 = face->vertex(1)->point();
    Point p2 = face->vertex(2)->point();

    K::FT maxEdge = CGAL::squared_distance(p0, p1);
    maxEdge = std::max(maxEdge, CGAL::squared_distance(p1, p2));
    maxEdge = std::max(maxEdge, CGAL::squared_distance(p2, p0));

    Point circumcenter = getCircumcenter(p0, p1, p2);
    K::FT circumradius = CGAL::squared_distance(circumcenter, p0);

    return circumradius / maxEdge;
}

void updatePheromones(std::map<CDT::Face_handle, double> &pheromones,
                      CDT &cdt,
                      const std::vector<Point> &cycleSteinerPoints,
                      double evaporationRate)
{
    // Evaporate pheromones globally
    for (auto &entry : pheromones)
    {
        entry.second *= (1.0 - evaporationRate);
        if (entry.second < 0.01)
        {
            entry.second = 0.01; // Avoid zero pheromones
        }
    }

    // Reinforce pheromones for improved faces
    for (const auto &steinerPoint : cycleSteinerPoints)
    {
        auto nearestFace = cdt.locate(steinerPoint); // Find the face containing the point
        if (nearestFace != cdt.infinite_face())
        {
            pheromones[nearestFace] += 1.0; // Reward this face
        }
    }
}

void antColonyOptimization(CDT &cdt, std::vector<Point> &steinerPoints,
                           int numAnts, int maxCycles, double evaporationRate)
{
    // Initialize pheromones for each finite face in CDT
    std::map<CDT::Face_handle, double> pheromones;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        pheromones[face] = 1.0; // Initial pheromone level for each face
    }

    for (int cycle = 0; cycle < maxCycles; ++cycle)
    {
        bool improvementMade = false;
        std::vector<Point> cycleSteinerPoints;

        for (int ant = 0; ant < numAnts; ++ant)
        {
            CDT tempCdt = cdt; // Work on a copy of the triangulation to evaluate improvements

            for (auto face = tempCdt.finite_faces_begin(); face != tempCdt.finite_faces_end(); ++face)
            {
                if (tempCdt.is_infinite(face))
                    continue; // Skip infinite faces

                // Check if the triangle formed by the face is obtuse
                if (isObtuse(face))
                {
                    // Generate candidate points to resolve the obtuse triangle
                    auto candidates = generateCandidatePoints(
                        face->vertex(0)->point(),
                        face->vertex(1)->point(),
                        face->vertex(2)->point());

                    // Try each candidate point for insertion
                    for (const auto &candidate : candidates)
                    {
                        tempCdt.insert(candidate); // Try inserting the candidate point into tempCdt

                        // Compare the number of obtuse triangles before and after inserting the candidate
                        int currentObtuse = countObtuseTriangles(tempCdt);
                        int globalObtuse = countObtuseTriangles(cdt);

                        if (currentObtuse < globalObtuse)
                        {
                            cdt = tempCdt;                           // Update the global CDT with the improved one
                            steinerPoints.push_back(candidate);      // Record the new Steiner point
                            cycleSteinerPoints.push_back(candidate); // Record this point for the current cycle
                            improvementMade = true;
                            break; // Break once an improvement is made (for this face)
                        }
                        else
                        {
                            // Find the nearest vertex manually (since `nearest_vertex` is not available)
                            auto nearestVertex = tempCdt.vertices_begin();
                            double minDistance = std::numeric_limits<double>::max();

                            // Iterate through vertices in tempCdt and find the nearest to the candidate
                            for (auto v = tempCdt.vertices_begin(); v != tempCdt.vertices_end(); ++v)
                            {
                                double dist = CGAL::squared_distance(candidate, v->point());
                                if (dist < minDistance)
                                {
                                    nearestVertex = v;
                                    minDistance = dist;
                                }
                            }

                            // Remove the nearest vertex if necessary
                            if (nearestVertex != tempCdt.vertices_end())
                            {
                                tempCdt.remove(nearestVertex);
                            }
                        }
                    }
                }
            }
        }

        if (cycleSteinerPoints.empty())
        {
            std::cout << "No Steiner points were inserted. Exiting loop." << std::endl;
            break; // Exit loop if no new Steiner points were inserted
        }

        std::cout << "Cycle " << cycle + 1
                  << ": Obtuse triangles remaining = " << countObtuseTriangles(cdt) << std::endl;

        // Update pheromones based on the improvements made during this cycle
        updatePheromones(pheromones, cdt, cycleSteinerPoints, evaporationRate);
    }
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
    int initialObtuseTriangles = countObtuseTriangles(cdt);
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
            int newObtuseTriangles = countObtuseTriangles(cdt);
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

std::vector<Point> performTriangulation(const json &inputData, CDT &cdt)
{
    // Validate input
    if (!inputData["points_x"].is_array() || !inputData["points_y"].is_array())
    {
        std::cerr << "Error: points_x and points_y must be arrays." << std::endl;
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

    // Step 2: Insert constraints into the CDT
    if (inputData.contains("additional_constraints"))
    {
        const auto &constraints = inputData["additional_constraints"];
        for (const auto &constraint : constraints)
        {
            size_t index1 = constraint[0];
            size_t index2 = constraint[1];

            Point p1(inputData["points_x"][index1], inputData["points_y"][index1]);
            Point p2(inputData["points_x"][index2], inputData["points_y"][index2]);
            cdt.insert_constraint(p1, p2);
        }
    }

    // Step 3: Choose triangulation method based on the "method" field
    std::string method = inputData.value("method", "local");
    const auto &params = inputData["parameters"];
    std::vector<Point> steinerPoints;

    if (method == "local")
    {
        int L = params.value("L", 1000);
        localSearchForSteinerPoints(cdt, steinerPoints, L);
    }
    else if (method == "simulated_annealing")
    {
        double alpha = params.value("alpha", 1.0);
        double beta = params.value("beta", 2.0);
        int L = params.value("L", 100);
        simulatedAnnealing(cdt, steinerPoints, alpha, beta, L);
    }
    else if (method == "aco")
    {
        double alpha = params.value("alpha", 1.0);
        double beta = params.value("beta", 2.0);
        double lambda = params.value("lambda", 0.7);
        int kappa = params.value("kappa", 10);
        int L = params.value("L", 100);
        antColonyOptimization(cdt, steinerPoints, kappa, L, lambda);
    }
    else
    {
        std::cerr << "Error: Unknown triangulation method '" << method << "'." << std::endl;
        exit(1);
    }

    // Extract triangulation results
    json result = extractTriangulationResults(cdt);
    result["num_points"] = cdt.number_of_vertices();
    result["num_edges"] = result["edges"].size();
    result["steiner_points_x"] = json::array();
    result["steiner_points_y"] = json::array();

    for (const auto &sp : steinerPoints)
    {
        result["steiner_points_x"].push_back(CGAL::to_double(sp.x()));
        result["steiner_points_y"].push_back(CGAL::to_double(sp.y()));
    }

    result["message"] = "Triangulation completed successfully!";
    result["method"] = method;

    return steinerPoints;
}

void writeOutput(const json &triangulationData, const CDT &cdt, const std::vector<Point> &steinerPoints, const std::string &filename)
{
    json outputData;

    // Copy input triangulation data
    outputData["content_type"] = "CG_SHOP_2025_Solution";
    outputData["instance_uid"] = triangulationData["instance_uid"];

    // Add Steiner points (if any) to the output data
    outputData["steiner_points_x"] = json::array();
    outputData["steiner_points_y"] = json::array();
    for (const auto &sp : steinerPoints)
    {
        outputData["steiner_points_x"].push_back(CGAL::to_double(sp.x()));
        outputData["steiner_points_y"].push_back(CGAL::to_double(sp.y()));
    }

    // Extract edges from the triangulation
    outputData["edges"] = json::array();
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face)
    {
        if (cdt.is_infinite(face))
            continue;

        // Add edges of the current triangle
        for (int i = 0; i < 3; ++i)
        {
            auto v1 = face->vertex(i);
            auto v2 = face->vertex((i + 1) % 3);

            // Find the index of each vertex
            int index1 = -1, index2 = -1;
            int idx = 0;
            for (auto v = cdt.vertices_begin(); v != cdt.vertices_end(); ++v)
            {
                if (v->point() == v1->point())
                    index1 = idx;
                if (v->point() == v2->point())
                    index2 = idx;
                ++idx;
            }

            if (index1 > index2)
                std::swap(index1, index2);

            outputData["edges"].push_back({index1, index2});
        }
    }

    // Add method and parameters (you can define them as needed)
    outputData["method"] = triangulationData["method"];
    outputData["parameters"] = triangulationData["parameters"];

    // Write the JSON to file
    std::ofstream outputFile(filename);
    if (!outputFile)
    {
        throw std::runtime_error("Could not open output file.");
    }
    outputFile << outputData.dump(4); // Pretty print the JSON with an indent of 4 spaces
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

void hybridOptimization(CDT &cdt, std::vector<Point> &steinerPoints, double alpha, double beta, int saIterations, int acoCycles, double evaporationRate)
{
    // Perform simulated annealing first
    simulatedAnnealing(cdt, steinerPoints, alpha, beta, saIterations);

    // Then perform ant colony optimization
    antColonyOptimization(cdt, steinerPoints, 10 /*numAnts*/, acoCycles, evaporationRate);
}

double evaluateMaxEdgeLength(const CDT &cdt)
{
    double maxEdge = 0.0;
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge)
    {
        auto v1 = edge->first->vertex(edge->second)->point();
        auto v2 = edge->first->vertex(edge->third)->point();
        maxEdge = std::max(maxEdge, CGAL::squared_distance(v1, v2));
    }
    return maxEdge;
}

double evaluateTriangulationWithEdgeLength(const CDT &cdt)
{
    double obtusePenalty = evaluateTriangulation(cdt); // Existing obtuse triangle score
    double maxEdgeLength = evaluateMaxEdgeLength(cdt);
    return obtusePenalty + 0.1 * maxEdgeLength; // Weighted combination
}

Point generateRandomPointInTriangle(const Point &p0, const Point &p1, const Point &p2)
{
    double r1 = ((double)rand() / RAND_MAX);
    double r2 = ((double)rand() / RAND_MAX);
    if (r1 + r2 > 1)
    {
        r1 = 1 - r1;
        r2 = 1 - r2;
    }
    double x = p0.x() + r1 * (p1.x() - p0.x()) + r2 * (p2.x() - p0.x());
    double y = p0.y() + r1 * (p1.y() - p0.y()) + r2 * (p2.y() - p0.y());
    return Point(x, y);
}
