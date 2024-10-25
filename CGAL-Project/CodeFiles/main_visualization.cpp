#include <QApplication>
#include "../HeaderFiles/visualization.h"
#include "../HeaderFiles/triangulation_utils.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // Create the CDT and define points
    CDT cdt;
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

    // Create visualization window
    VisualizationWindow window;
    window.displayTriangulation(cdt);
    window.resize(800, 600);
    window.show();

    // Insert Steiner points and apply diagonal flips to eliminate obtuse triangles
    insertSteinerPoints(cdt);

    // Display refined triangulation
    window.displayTriangulation(cdt);

    return app.exec();
}
