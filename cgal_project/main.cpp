
// Mock main program to showcase the the usage of the functions
// To run 
// -cd built 
// -cmake ..
// -make 
// -./triangulation_program

#include "triangulation_utils.h"
#include <iostream>

int main() {
    CDT cdt;

    // Define points that form an obtuse triangle
    std::vector<Point> points = {
        Point(0, 0),   // A
        Point(4, 0),   // B
        Point(1, 2),   // C (angle at this point is obtuse)
        Point(3, 2)
    };

    // Insert initial points into CDT
    for (const auto& p : points) {
        cdt.insert(p);
    }

    std::cout << "Initial triangulation:\n";
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        auto v1 = edge->first->vertex((edge->second + 1) % 3);
        auto v2 = edge->first->vertex((edge->second + 2) % 3);
        std::cout << "Edge (" << v1->point() << ") - (" << v2->point() << ")\n";
    }

    // Insert Steiner points to eliminate obtuse triangles
    insertSteinerPoints(cdt);

    std::cout << "\nTriangulation after inserting Steiner points:\n";
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        auto v1 = edge->first->vertex((edge->second + 1) % 3);
        auto v2 = edge->first->vertex((edge->second + 2) % 3);
        std::cout << "Edge (" << v1->point() << ") - (" << v2->point() << ")\n";
    }

    return 0;
}
