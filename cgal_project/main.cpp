#include "triangulation_utils.h"
#include <iostream>

int main()
{
    CDT cdt;

    // Define points that form an obtuse triangle and additional points
    std::vector<Point> points = {
        Point(0, 0), // A
        Point(1, 0), // B
        Point(5, 1), // C (angle at this point is obtuse)
        Point(8, 6), // D
        Point(2, 6), // E
        Point(1, 5), // F
        Point(3, 4), // G
        Point(8, 8)  // H
    };

    // Insert initial points into CDT
    for (const auto &p : points)
    {
        cdt.insert(p);
    }

    std::cout << "Initial triangulation:\n";
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge)
    {
        auto v1 = edge->first->vertex((edge->second + 1) % 3);
        auto v2 = edge->first->vertex((edge->second + 2) % 3);
        std::cout << "Edge (" << v1->point() << ") - (" << v2->point() << ")\n";
    }

    // Insert Steiner points and apply diagonal flips to eliminate obtuse triangles
    std::cout << "\nRefining triangulation to remove obtuse angles...\n";
    insertSteinerPoints(cdt);

    std::cout << "\nTriangulation after inserting Steiner points and applying diagonal flips:\n";
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge)
    {
        auto v1 = edge->first->vertex((edge->second + 1) % 3);
        auto v2 = edge->first->vertex((edge->second + 2) % 3);
        std::cout << "Edge (" << v1->point() << ") - (" << v2->point() << ")\n";
    }

    return 0;
}