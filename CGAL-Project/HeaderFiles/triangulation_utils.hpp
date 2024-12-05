#ifndef TRIANGULATION_UTILS_HPP
#define TRIANGULATION_UTILS_HPP

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

// Function prototypes
void parseInput(const string &inputFile, json &inputData);
bool isObtuse(const CDT::Face_handle &face);
int countObtuseTriangles(CDT &cdt, const Point &steinerPoint);
void localSearchForSteinerPoints(CDT &cdt, vector<Point> &steinerPoints, int maxIterations);
json performTriangulation(const json &inputData, CDT &cdt);
void writeOutput(const json &inputData, const CDT &cdt, const vector<Point> &steinerPoints, const string &outputFile);
json extractTriangulationResults(const CDT &cdt);

#endif // TRIANGULATION_UTILS_HPP
