#ifndef TRIANGULATION_UTILS_HPP
#define TRIANGULATION_UTILS_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <vector>
#include <nlohmann/json.hpp>

// Kernel and types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point;
typedef K::Vector_2 Vector;
typedef K::FT FT;
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Face_handle Face_handle;

// Edge type definition (pair of face handle and index)
typedef std::pair<Face_handle, int> Edge;

// Function declarations
void parseInput(const nlohmann::json &input, std::vector<Point> &points, std::vector<std::pair<int, int>> &constraints);
bool isObtuse(CDT::Face_handle face);
bool isObtuseAngle(const Point &p1, const Point &p2, const Point &p3);
Face_handle findObtuseAdjacentFace(Face_handle face);
void writeOutput(const CDT &cdt, const std::vector<Point> &steinerPoints, nlohmann::json &output);
bool applyEdgeFlip(CDT &cdt, CDT::Face_handle &face);
void insertSteinerPoint(CDT &cdt, CDT::Face_handle face);
Point insertSteinerForPolygon(const std::vector<Point> &polygon_points);
void applyPolygonalTreatment(CDT &cdt, CDT::Face_handle face1, CDT::Face_handle face2);
void flipEdgeIfBeneficial(CDT &cdt, const Edge &edge); // Add function declaration here

#endif // TRIANGULATION_UTILS_HPP
