#ifndef TRIANGULATION_UTILS_H
#define TRIANGULATION_UTILS_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds, Itag> CDT;
typedef CDT::Point Point;

bool isObtuse(const CDT::Face_handle& face);

Point calculateSteinerPoint(const CDT& cdt);
void insertSteinerPoints(CDT& cdt);
bool tryDiagonalFlip(CDT& cdt , CDT::Face_handle& face);

#endif // TRIANGULATION_UTILS_H
